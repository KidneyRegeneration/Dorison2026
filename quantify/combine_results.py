import pandas as pd
import argparse
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
                       help='input sample names, if not specified, csv filenames will be used.')
    
    parser.add_argument('-o', '--output_csv',
                       required=True, 
                       type=str,
                       help='output csv filepath')
    
    parser.add_argument('-n', '--stain_name',
                       required=False, 
                       type=str,
                       default="KAT2B",
                       help='Name of stain')

    
    return parser.parse_args()

def combine_csvs(csv_list, output_csv, samplenames=None, stain_name="KAT2B"):
    mean_values = {}
    for i, csv in enumerate(csv_list):
        df = pd.read_csv(csv)
        mean_df = pd.DataFrame(df.mean()).T.reset_index(drop=True)
        total_df = pd.DataFrame(df.sum()).T.reset_index(drop=True)
        merged_df = pd.concat([mean_df.add_suffix('_avg'), total_df.add_suffix('_sum')], axis=1)

        if samplenames is not None:
            sample = samplenames[i]
        else:
            sample = Path(csv).stem
        metadata = sample.split("_")
        metadata[-2] = metadata[-2] + "_" + metadata[-1]
        metadata = metadata[:-1]
        for metadata_name, val in zip(["stain", "condition", "experiment", "replicate", "if_exp", "unique_id"], metadata):
            merged_df[metadata_name] = val
        merged_df.drop(labels=["dapi_min_sum", "dapi_max_sum", "dapi_mean_sum", "stain_mean_in_pods_sum", "stain_mean_out_pods_sum", "stain_mean_out_nuclei_sum"], axis=1, inplace=True)
        merged_df.rename({"dapi_mean_avg": "dapi_mean",
                          "dapi_min_avg": "dapi_min",
                          "dapi_max_avg": "dapi_max",
                          "stain_mean_in_pods_avg": f"{stain_name}_mean_in_pods",
                          "stain_mean_out_pods_avg": f"{stain_name}_mean_out_pods",
                          "stain_mean_out_nuclei_avg": f"{stain_name}_mean_out_nuclei",
                          }, axis=1, inplace=True)

        merged_df[f"{stain_name}_per_pod"] = merged_df["stain_sum_in_pods_sum"] / merged_df["pod_size_sum"]
        merged_df[f"{stain_name}_sum_norm_size_avg_div_dapi_mean"] = merged_df["stain_sum_norm_size_avg"] / merged_df["dapi_mean"]
        merged_df[f"{stain_name}_in_out_pods"] = merged_df[f"{stain_name}_mean_in_pods"] - merged_df[f"{stain_name}_mean_out_pods"]
        merged_df[f"{stain_name}_in_out_nuclei"] = merged_df[f"{stain_name}_mean_in_pods"] - merged_df[f"{stain_name}_mean_out_nuclei"]
        mean_values[sample] = merged_df
    
    df = pd.concat(mean_values)
    df.columns = df.columns.str.replace('stain', stain_name, case=False)
    # df = pd.DataFrame([mean_values]).T
    df = df.dropna()
    df.to_csv(output_csv)

def main():
    args = parse_args()
    combine_csvs(args.input_csvs, args.output_csv, samplenames=args.input_samplenames, stain_name=args.stain_name)

if __name__ == "__main__":
    main()