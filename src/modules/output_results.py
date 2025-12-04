import matplotlib.pyplot as plt
import numpy as np

def save_plot(strand, predictions, window_size, step, figure_folder):
    """ Plots the genomic sequence. """
    plt.figure(figsize=(14, 6))

    plt.rcParams.update({
        'font.size': 16,
        'axes.titlesize': 22,
        'axes.labelsize': 18,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14,
    })
    
    positions = [(start + end) / 2 for (start, end) in sorted(predictions.keys())]
    probabilities = [predictions[(start, end)] for (start, end) in sorted(predictions.keys())]
    
    plt.plot(positions, probabilities, label="P(genic)", linewidth=2)

    plt.title(f'Strand: {strand} - Predictions', pad=15)
    plt.xlabel('Genomic Position (center of window)')
    plt.ylabel('Probability of being genic')
    plt.legend()
    plt.grid(True, linewidth=0.6, alpha=0.6)
    plt.tight_layout()
    
    save_path = f"{figure_folder}/{strand}_combined_{window_size}_{step}.png"
    plt.savefig(save_path)
    print(f"[INFO] Combined plot saved to {save_path}")
    plt.close()

def merge_overlapping(intervals, scores):
    """ Merges overlapping intervals and averages their scores. """

    if not intervals: 
        return []

    merged = []
    cs, ce, cscores = intervals[0][0], intervals[0][1], [scores[0]]

    for (s, e), sc in zip(intervals[1:], scores[1:]):
        if s <= ce:
            ce = max(ce, e)
            cscores.append(sc)
        else:
            merged.append((cs, ce, np.mean(cscores)))
            cs, ce, cscores = s, e, [sc]
    merged.append((cs, ce, np.mean(cscores)))

    return merged

def print_predictions(predictions, chromosome, strand, threshold, out_file, \
                        source="GeAnno", append=False, start_id=1):
    """ Prints predicted genic regions to the output file. """

    chromosome = str(chromosome).lstrip(">")
    intervals, scores, in_block = [], [], False

    for (s, e) in sorted(predictions):
        prob = predictions[(s, e)]
        if prob >= threshold:
            if not in_block:
                block_start = s
            block_end, last_prob, in_block = e, prob, True
        else:
            if in_block:
                intervals.append((block_start, block_end))
                scores.append(last_prob)
                in_block = False

    if in_block:
        intervals.append((block_start, block_end))
        scores.append(last_prob)

    merged = merge_overlapping(intervals, scores)

    gff_strand = '+' if (strand == 'forward' or strand == "1") else '-' if (strand == 'reverse' or strand == "0") else strand

    written = 0
    i = start_id
    mode = "a" if append else "w"
    with open(out_file, mode) as f:
        if not append:
            f.write("##gff-version 3\n")
        
        for s, e, m in merged:
            attributes = f"ID=gene{i};"
            f.write(f"{chromosome}\t{source}\tgene\t{s}\t{e}\t{m:.4f}\t{gff_strand}\t.\t{attributes}\n")
            i += 1
            written += 1

    return written