import os, logging, joblib, argparse
import numpy as np
import pandas as pd
from collections import defaultdict

from modules.output_results import save_plot, print_predictions
from modules.dna_import import reverse_complement, import_dna_file
from modules.subprocesses import run_orfs, run_jarvis, get_all_features

logging.basicConfig(level=logging.INFO, format='[PRED] %(message)s')

MOTIFS = "A,C,T,G,ATG,TAA,TAG,TGA,TATAAA,CAAT,CGCG,GT,AG,AAA,TTT,GGG,CGG,CAG"

class PipelineRun:
    def __init__(self, dna, model, window_size, step, thresholds, moving_average,
                 csv_input, output_folder, figure_folder, features_folder):
        
        thresholds_float = [float(x) for x in thresholds.split(",")]

        if figure_folder is not None:
            os.makedirs(figure_folder, exist_ok=True)
        os.makedirs(output_folder, exist_ok=True)
        
        features_file = f"all_windows_{window_size}_{step}.csv"

        if not features_folder:
            features_path = None
        else:
            os.makedirs(features_folder, exist_ok=True)
            features_path = os.path.join(features_folder, features_file)

        self.dna = dna
        self.model = model
        self.window_size = window_size
        self.step = step
        self.thresholds = thresholds_float
        self.moving_average = moving_average
        self.csv_input = csv_input
        self.output_folder = output_folder
        self.figure_folder = figure_folder
        self.features_path = features_path


def prepare_windows(chr_dict, window_size, step):
    """ Precompute all windows first, without ORFs or compressions.
        Returns list of windows and a dict mapping (chromosome, strand) to list of (start, end) tuples.
    """

    windows_by_chrom_strand = defaultdict(list)
    windows = []

    for chromosome, chr_seq in chr_dict.items():
        for strand in ("forward", "reverse"):
            seq = chr_seq if strand == "forward" else reverse_complement(chr_seq)
            seq_len = len(seq)
            for start in range(0, len(seq) - window_size + 1, step):
                end = start + window_size
                seq_window = seq[start:end]
                windows.append((chromosome, strand, start, end, seq_window))
                windows_by_chrom_strand[(chromosome, strand)].append((start, end))

            # Handle last window if it doesn't fit perfectly
            if end < seq_len:
                start = seq_len - window_size
                seq_window = seq[start:seq_len]
                windows.append((chromosome, strand, start, seq_len, seq_window))
                windows_by_chrom_strand[(chromosome, strand)].append((start, seq_len))

    return windows, windows_by_chrom_strand


def apply_moving_average(predictions, ma_window_size=5):
    """ Applies a moving average to the predicted probabilities to smooth the results """

    sorted_windows = sorted(predictions.keys(), key=lambda x: (x[0] + x[1]) / 2)
    scores = [predictions[window] for window in sorted_windows]
    smoothed_scores = []
    half_window = ma_window_size // 2

    for i in range(len(scores)):
        start_idx = max(0, i - half_window)
        end_idx = min(len(scores), i + half_window + 1)
        smoothed_scores.append(np.mean(scores[start_idx:end_idx]))

    return dict(zip(sorted_windows, smoothed_scores))


def extract_windows_and_features(chr_dict, pipeline):
    """ Generate windows, compute ORFs and compressions, and extract all features.
        Returns the final DataFrame and a dict mapping (chromosome, strand) to list of (start, end) tuples.
    """

    windows, windows_by_chrom_strand = prepare_windows(chr_dict, pipeline.window_size, pipeline.step)

    orfs = run_orfs(windows)
    compressions = run_jarvis(windows)

    final_df = get_all_features(windows, orfs, compressions, MOTIFS, pipeline.features_path)
    
    return final_df, windows_by_chrom_strand


def run_pipeline(pipeline: PipelineRun):
    """ Main function to run the genic prediction pipeline """
    
    model_job = joblib.load(pipeline.model)
    model = model_job["pipeline"]
    le = model_job["label_encoder"]
    feature_names = model_job["feature_names"]

    if pipeline.csv_input and os.path.exists(pipeline.csv_input):
        logging.info(f"[MAIN] Loading pre-computed features from {pipeline.csv_input}.")
        final_df = pd.read_csv(pipeline.csv_input)
        final_df["strand"] = final_df["strand"].map({1: "forward", 0: "reverse"}).fillna("forward")

        windows_by_chrom_strand = defaultdict(list)
        for _, row in final_df.iterrows():
            row["strand"]
            windows_by_chrom_strand[(row["chromosome"], row["strand"])].append((row["start"], row["end"]))

        chr_lengths = final_df.groupby("chromosome")["end"].max().to_dict()

    else:
        logging.info("[MAIN] Generating windows and calculating characteristics.")
        chr_dict = import_dna_file(pipeline.dna)
        final_df, windows_by_chrom_strand = extract_windows_and_features(chr_dict, pipeline)
        chr_lengths = {chrom: len(seq) for chrom, seq in chr_dict.items()}

    features_df = final_df.loc[:, ~final_df.columns.duplicated()]
    features_df = features_df[[col for col in feature_names if col in features_df.columns]]

    logging.info(f"[MAIN] Final features extracted.")

    p_genic = model.predict_proba(features_df)[:, list(le.classes_).index("gene")]

    prediction_cache = defaultdict(dict)
    i = 0
    append=False
    for (chromosome, strand), windows_list in windows_by_chrom_strand.items():
        chr_len = chr_lengths[chromosome]
        for start, end in windows_list:
            if (strand == "reverse"):
                start_m, end_m = chr_len - end, chr_len - start
            else:
                start_m, end_m = start, end
            prediction_cache[(chromosome, strand)][(start_m, end_m)] = p_genic[i]
            i += 1

    initialized_files = set()
    next_gene_id = defaultdict(lambda: 1)
    for (chromosome, strand), predictions in prediction_cache.items():
        smoothed = apply_moving_average(predictions, pipeline.moving_average)

        if pipeline.figure_folder is not None:
            save_plot(strand, smoothed, pipeline.window_size, pipeline.step, pipeline.figure_folder)
        
        for threshold in pipeline.thresholds:
            out_file = os.path.join(
                pipeline.output_folder, 
                f"output_{pipeline.window_size}_{pipeline.step}_{threshold}.gff3"
            )

            append = out_file in initialized_files
            written = print_predictions(smoothed, chromosome, strand, threshold, out_file, source="GeAnno", append=append, start_id=next_gene_id[out_file])

            initialized_files.add(out_file)
            next_gene_id[out_file] += written

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    #required parameters
    parser.add_argument("-d", "--dna_file", help="FASTA DNA file (required if --csv_features_file_path not provided)")
    parser.add_argument("-m", "--model", required=True)
    
    #where to place outputs
    parser.add_argument("-f", "--figure_folder", default="./genic/figures")
    parser.add_argument("-o", "--output_folder", default="./genic/output")
    parser.add_argument("-ff", "--features_folder", default=None)

    #where to search for already computed features
    parser.add_argument("-c", "--csv_features_file_path", default="") 
 
    #optional parameters
    parser.add_argument("-p", "--save_plot", type=bool, default=False)
    parser.add_argument("-s", "--step", default=50, type=int)
    parser.add_argument("-t", "--thresholds", default="0.8", type=str)
    parser.add_argument("-w", "--window_size", default=1500, type=int)
    parser.add_argument("-ma", "--moving_average", default=5, type=int)
    
    args = parser.parse_args()

    if not args.csv_features_file_path and not args.dna_file:
        parser.error("Either --dna_file (-d) or --csv_features_file_path (-c) must be provided.")

    if not (args.save_plot):
        args.figure_folder = None
        print("[INFO] Plots will not be saved as per user request. Ignoring --figure_folder argument.")

    pipeline = PipelineRun(args.dna_file, args.model, args.window_size, args.step, 
                 args.thresholds, args.moving_average, args.csv_features_file_path,
                 args.output_folder, args.figure_folder, args.features_folder)

    run_pipeline(pipeline)