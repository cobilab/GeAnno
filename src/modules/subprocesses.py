import os, re, subprocess, tempfile, logging
import pandas as pd
from io import StringIO
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_orfs(windows):
    """ Run getorf on each window in parallel and return list of longest ORFs. """
    
    orfs = [None] * len(windows)

    def orf_job(idx, seq):
        try:
            with tempfile.NamedTemporaryFile(delete=False, suffix=".seq", mode="w") as f:
                f.write(seq)
                seq_file = f.name
            orf_file = seq_file.replace(".seq", ".orf")
            
            subprocess.run(
                ["getorf", "-sequence", seq_file, "-outseq", orf_file, "-reverse", "no"],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )

            longest, current = "", []
            with open(orf_file) as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if line.startswith(">"):
                        if len(current) > len(longest):
                            longest = "".join(current)
                        current = []
                    else:
                        current.append(line)
                if len(current) > len(longest):
                    longest = "".join(current)

        except subprocess.CalledProcessError:
            longest = ""
        finally:
            for fname in (seq_file, orf_file):
                if os.path.exists(fname):
                    os.remove(fname)

        return idx, longest

    logging.info("[MAIN] Running ORF extraction in parallel")
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(orf_job, i, win[4]) for i, win in enumerate(windows)]
        for future in tqdm(as_completed(futures), total=len(futures), desc="ORFs", unit="win"):
            idx, orf = future.result()
            orfs[idx] = orf

    return orfs

def run_jarvis(windows):
    """ Run JARVIS3 on each window in parallel and return list of compression ratios. """

    compressions = [None] * len(windows)

    def jarvis_job(idx, seq):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".seq", mode="w", dir=".") as f:
            f.write(seq)
            seq_filename = f.name

        result = subprocess.run(
            ['JARVIS3', '-t', '4', '-lr', '0', '-cm', '1:1:0:0.9/0:0:0:0', seq_filename],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )

        nc = "0.0"
        if result.returncode == 0:
            match = re.search(r'([0-9.]+) NC', result.stdout)
            if match:
                nc = match.group(1)
        else:
            print(f"=== JARVIS3 ERROR on idx {idx} ===")
            print(f"stderr: {result.stderr}")
            print(f"stdout: {result.stdout}")
            print(f"sequence file: {seq_filename}")

        os.remove(seq_filename)
        jc_filename = seq_filename + ".jc"
        if os.path.exists(jc_filename):
            os.remove(jc_filename)

        return (idx, nc)

    logging.info("[MAIN] Running JARVIS3 in parallel")
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(jarvis_job, i, win[4]) for i, win in enumerate(windows)]
        for future in tqdm(as_completed(futures), total=len(futures), desc="JARVIS3 windows", unit="win"):
            idx, comp = future.result()
            compressions[idx] = comp

    return compressions

def get_all_features(windows, orfs, compressions, motifs, output_csv=None):
    """
    Write temp .tsv file from sequences, ORFs, and compressions.
    Run external binary and return full feature DataFrame including metadata.
    """
    cpp_input_df = pd.DataFrame([
        [seq, orfs[i], compressions[i], 1 if strand == "forward" else 0]
        for i, (_, strand, _, _, seq) in enumerate(windows)
    ], columns=["sequence", "orf", "compression", "strand"])

    with tempfile.NamedTemporaryFile(delete=False, suffix=".tsv", mode="w") as f:
        tmp_tsv = f.name
        cpp_input_df.to_csv(f, sep="\t", index=False, header=False)

    result = subprocess.run(['extract_characteristics_batch', '-i', tmp_tsv, '-t', motifs],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
    os.remove(tmp_tsv)

    features_df = pd.read_csv(StringIO(result.stdout))

    meta_df = pd.DataFrame([
        [chrom, start, end, seq]
        for (chrom, _, start, end, seq) in windows
    ], columns=["chromosome", "start", "end", "sequence"])

    final_df = pd.concat([meta_df, features_df], axis=1)

    if output_csv:
        # save to csv so that it can be saved later
        final_df.to_csv(output_csv, index=False)

    return features_df
