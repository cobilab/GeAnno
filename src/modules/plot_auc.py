import os
import pandas as pd
import matplotlib.pyplot as plt

def plot_curves_in_folder(folder):
    files = os.listdir(folder)
    
    roc_files = [f for f in files if f.endswith('_roc.csv')]
    
    for roc_file in roc_files:
        prefix = roc_file[:-8]  # _roc.csv
        prc_file = prefix + '_prc.csv'
        
        roc_path = os.path.join(folder, roc_file)
        prc_path = os.path.join(folder, prc_file)
        
        if os.path.exists(roc_path):
            df_roc = pd.read_csv(roc_path)
            plt.figure()
            plt.plot(df_roc['FPR'], df_roc['TPR'], label='ROC curve')
            plt.plot([0, 1], [0, 1], 'k--', label='Random')
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title(f'ROC Curve: {prefix}')
            plt.legend()
            plt.grid(True)
            plt.savefig(os.path.join(folder, f'{prefix}_roc.png'))
            plt.close()
        
        if os.path.exists(prc_path):
            df_prc = pd.read_csv(prc_path)
            plt.figure()
            plt.plot(df_prc['Recall'], df_prc['Precision'], label='PRC curve')
            plt.xlabel('Recall')
            plt.ylabel('Precision')
            plt.title(f'Precision-Recall Curve: {prefix}')
            plt.legend()
            plt.grid(True)
            plt.savefig(os.path.join(folder, f'{prefix}_prc.png'))
            plt.close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Plot ROC and PRC curves from CSV files in a folder.")
    parser.add_argument('folder', help="Folder containing the CSV files")
    args = parser.parse_args()
    
    plot_curves_in_folder(args.folder)
