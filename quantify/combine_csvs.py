import pandas as pd
import argparse
import json
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Get stain amounts in each podocyte')
    parser.add_argument('-i', '--input_csvs',
                        required=True,
                        nargs="+",
                        help='Input csv filepaths')
    parser.add_argument('-s', '--input_samplenames',
                        required=False,
                        default=None,
                        nargs="+",
                        help='input sample names, if not specified, csv filenames will be used.')
    parser.add_argument('-m', '--input_metadata',
                        required=False,
                        default=None,
                        nargs="+",
                        help='input metadata json files, one per csv. If provided, metadata is read from these files.')
    parser.add_argument('-o', '--output_csv',
                        required=True,
                        type=str,
                        help='output csv filepath')
    return parser.parse_args()

def combine_csvs(csv_list, output_csv, samplenames=None, metadata_files=None):
    mean_values = {}
    for i, csv in enumerate(csv_list):
        df = pd.read_csv(csv)
        mean_df = pd.DataFrame(df.mean(numeric_only=True)).T.reset_index(drop=True)

        if metadata_files is not None:
            # Read metadata from JSON file
            with open(metadata_files[i]) as f:
                meta = json.load(f)
            sample = meta.get('sample_name', Path(csv).stem)
            mean_df['stain']      = meta.get('stain', '')
            mean_df['condition']  = meta.get('condition', '')
            mean_df['condition2']  = meta.get('condition2', '')
            mean_df['experiment'] = meta.get('experiment', '')
            mean_df['replicate']  = meta.get('replicate', '')
            mean_df['if_exp']     = meta.get('if_exp', '')
            mean_df['unique_id']  = meta.get('unique_id', '')
        else:
            # Fall back to sample name parsing
            if samplenames is not None:
                sample = samplenames[i]
            else:
                sample = Path(csv).stem
            metadata = sample.split("_")
            metadata[-2] = metadata[-2] + "_" + metadata[-1]
            metadata = metadata[:-1]
            for metadata_name, val in zip(["stain", "condition", "experiment", "replicate", "if_exp", "unique_id"], metadata):
                mean_df[metadata_name] = val

        mean_values[sample] = mean_df

    df = pd.concat(mean_values)
    df = df.dropna()
    if "podocyte_index" in df.columns:
        df.drop(labels=["podocyte_index"], axis=1, inplace=True)
    df.to_csv(output_csv)

def main():
    args = parse_args()
    combine_csvs(args.input_csvs, args.output_csv,
                 samplenames=args.input_samplenames,
                 metadata_files=args.input_metadata)

if __name__ == "__main__":
    main()