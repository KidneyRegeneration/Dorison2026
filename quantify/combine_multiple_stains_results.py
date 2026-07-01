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
    
    return parser.parse_args()


def combine_csvs(csv_list, output_csv, samplenames=None):
    mean_values = {}
    
    for i, csv in enumerate(csv_list):
        df = pd.read_csv(csv)
        
        # Detect stain name(s) from columns by looking for known suffixes
        stain_name = None
        known_suffixes = ['_mean_in_pods', '_mean_out_pods', '_mean_out_nuclei', 
                          '_sum_in_pods', '_sum_norm_size']
        
        for col in df.columns:
            for suffix in known_suffixes:
                if col.endswith(suffix):
                    stain_name = col.replace(suffix, '')
                    break
            if stain_name:
                break
        
        # Fallback to 'stain' if no match found
        if stain_name is None:
            stain_name = 'stain'
        
        mean_df = pd.DataFrame(df.mean()).T.reset_index(drop=True)
        total_df = pd.DataFrame(df.sum()).T.reset_index(drop=True)
        merged_df = pd.concat([mean_df.add_suffix('_avg'), total_df.add_suffix('_sum')], axis=1)
        
        if samplenames is not None:
            sample = samplenames[i]
        else:
            sample = Path(csv).stem
            
        metadata = sample.split("_")
        print("Metadata = ", metadata)
        metadata[-2] = metadata[-2] + "_" + metadata[-1]
        print("Metadata = ", metadata)
        metadata = metadata[1:-1]
        print("Metadata = ", metadata)
        for metadata_name, val in zip(["condition", "experiment", "replicate", "unique_id"], metadata):
            merged_df[metadata_name] = val
        
        # Drop columns with dynamic stain name detection
        cols_to_drop = [
            "dapi_min_sum", "dapi_max_sum", "dapi_mean_sum",
            f"{stain_name}_mean_in_pods_sum",
            f"{stain_name}_mean_out_pods_sum",
            f"{stain_name}_mean_out_nuclei_sum",
            f"podocyte_index_avg"
        ]

        merged_df.drop(labels=[col for col in cols_to_drop if col in merged_df.columns], 
                       axis=1, inplace=True)
        
        # Rename columns dynamically
        rename_dict = {
            "dapi_mean_avg": "dapi_mean",
            "dapi_min_avg": "dapi_min",
            "dapi_max_avg": "dapi_max",
            f"{stain_name}_mean_in_pods_avg": f"{stain_name}_mean_in_pods",
            f"{stain_name}_mean_out_pods_avg": f"{stain_name}_mean_out_pods",
            f"{stain_name}_mean_out_nuclei_avg": f"{stain_name}_mean_out_nuclei",
        }
        merged_df.rename({k: v for k, v in rename_dict.items() if k in merged_df.columns}, 
                         axis=1, inplace=True)
        
        # Calculate derived columns
        if f"{stain_name}_sum_in_pods_sum" in merged_df.columns and "pod_size_sum" in merged_df.columns:
            merged_df[f"{stain_name}_per_pod"] = merged_df[f"{stain_name}_sum_in_pods_sum"] / merged_df["pod_size_sum"]
        
        if f"{stain_name}_sum_norm_size_avg" in merged_df.columns and "dapi_mean" in merged_df.columns:
            merged_df[f"{stain_name}_sum_norm_size_avg_div_dapi_mean"] = merged_df[f"{stain_name}_sum_norm_size_avg"] / merged_df["dapi_mean"]
        
        if f"{stain_name}_mean_in_pods" in merged_df.columns and f"{stain_name}_mean_out_pods" in merged_df.columns:
            merged_df[f"{stain_name}_in_out_pods"] = merged_df[f"{stain_name}_mean_in_pods"] - merged_df[f"{stain_name}_mean_out_pods"]
        
        if f"{stain_name}_mean_in_pods" in merged_df.columns and f"{stain_name}_mean_out_nuclei" in merged_df.columns:
            merged_df[f"{stain_name}_in_out_nuclei"] = merged_df[f"{stain_name}_mean_in_pods"] - merged_df[f"{stain_name}_mean_out_nuclei"]
        
        mean_values[sample] = merged_df
    
    df = pd.concat(mean_values)
    df = df.dropna()
    df.to_csv(output_csv)


def main():
    args = parse_args()
    combine_csvs(args.input_csvs, args.output_csv, samplenames=args.input_samplenames)

if __name__ == "__main__":
    main()