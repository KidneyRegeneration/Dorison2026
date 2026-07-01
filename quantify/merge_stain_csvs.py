#!/usr/bin/env python3
import argparse
import pandas as pd

import pandas as pd
import numpy as np

def merge_with_mean(df1, df2, key_column='sample_name'):
    """
    Merge two dataframes on a key column, keeping all columns from both.
    For columns that exist in both dataframes (except the key), take the mean of values.
    
    Parameters:
    -----------
    df1 : pd.DataFrame
        First dataframe
    df2 : pd.DataFrame
        Second dataframe
    key_column : str
        Column name to merge on (default: 'sample_name')
    
    Returns:
    --------
    pd.DataFrame
        Merged dataframe with averaged values for overlapping columns
    """
    # Perform outer merge to keep all rows
    merged = pd.merge(df1, df2, on=key_column, how='outer', suffixes=('_df1', '_df2'))
    
    # Find columns that exist in both dataframes (excluding the key column)
    cols_df1 = set(df1.columns) - {key_column}
    cols_df2 = set(df2.columns) - {key_column}
    common_cols = cols_df1.intersection(cols_df2)
    
    # Process common columns: average values and remove suffixed versions
    for col in common_cols:
        col_df1 = f"{col}_df1"
        col_df2 = f"{col}_df2"
        
        
        # Check if columns are numeric
        is_numeric_df1 = pd.api.types.is_numeric_dtype(merged[col_df1])
        is_numeric_df2 = pd.api.types.is_numeric_dtype(merged[col_df2])
        
        if is_numeric_df1 and is_numeric_df2:
            # Calculate mean for numeric columns, handling NaN values appropriately
            merged[col] = merged[[col_df1, col_df2]].mean(axis=1)
        else:
            # Handle non-numeric columns
            merged[col] = None  # Initialize column
            for idx in merged.index:
                val1 = merged.loc[idx, col_df1]
                val2 = merged.loc[idx, col_df2]
                
                # Both are NaN/None - keep as NaN
                if pd.isna(val1) and pd.isna(val2):
                    merged.loc[idx, col] = np.nan
                # Only one value exists - use it
                elif pd.isna(val1):
                    merged.loc[idx, col] = val2
                elif pd.isna(val2):
                    merged.loc[idx, col] = val1
                # Both exist - check if they match
                elif val1 == val2:
                    merged.loc[idx, col] = val1  # Take first instance
                else:
                    # Values don't match - print error
                    print(f"ERROR: Mismatched non-numeric values in column '{col}' at sample '{merged.loc[idx, key_column]}': '{val1}' vs '{val2}'")
                    merged.loc[idx, col] = val1  # Take first instance despite mismatch
        
        # Drop the suffixed columns
        merged = merged.drop(columns=[col_df1, col_df2])
    
    # Process columns unique to each dataframe
    # These will have suffixes that need to be removed
    unique_df1_cols = cols_df1 - common_cols
    unique_df2_cols = cols_df2 - common_cols
    
    rename_dict = {}
    for col in unique_df1_cols:
        suffixed_name = f"{col}_df1"
        if suffixed_name in merged.columns:
            rename_dict[suffixed_name] = col
    
    for col in unique_df2_cols:
        suffixed_name = f"{col}_df2"
        if suffixed_name in merged.columns:
            rename_dict[suffixed_name] = col
    
    merged = merged.rename(columns=rename_dict)
    
    return merged


def main():
    parser = argparse.ArgumentParser(description='Merge two stain quantification CSV files')
    parser.add_argument('--stainA_csv', required=True, help='Path to stain A CSV file')
    parser.add_argument('--stainB_csv', required=True, help='Path to stain B CSV file')
    parser.add_argument('--stainA_name', required=True, help='Name of stain A')
    parser.add_argument('--stainB_name', required=True, help='Name of stain B')
    parser.add_argument('--output', required=True, help='Output merged CSV file path')
    parser.add_argument('--extra_input', required=False, help='Extra input for nextflow rerun')
    
    args = parser.parse_args()
    
    # Read the two CSV files
    df_a = pd.read_csv(args.stainA_csv)
    df_b = pd.read_csv(args.stainB_csv)
    formatted_string = ["_".join(i.split("_")[1:]) for i in df_a["Unnamed: 0"]]
    formatted_string_b = ["_".join(i.split("_")[1:]) for i in df_b["Unnamed: 0"]]
    
    df_a["sample_name"] = formatted_string
    df_b["sample_name"] = formatted_string_b
    
    df_a.drop("Unnamed: 0", axis=1, inplace=True)
    df_b.drop("Unnamed: 0", axis=1, inplace=True)
    
    # Replace "stain" in column names with the actual stain names
    df_a.columns = df_a.columns.str.replace('stain', args.stainA_name, case=False)
    df_b.columns = df_b.columns.str.replace('stain', args.stainB_name, case=False)
    
    df_merged = merge_with_mean(df_a, df_b)
    df_merged.index = df_merged["sample_name"]
    df_merged.to_csv(args.output)
    print(df_merged)
    exit()    
    # Find common columns (those that don't contain the stain names)
    cols_a_without_stain = [col for col in df_a.columns if args.stainA_name.lower() not in col.lower()]
    cols_b_without_stain = [col for col in df_b.columns if args.stainB_name.lower() not in col.lower()]
    
    # Get the intersection of these columns to use as merge keys
    merge_on = list(set(cols_a_without_stain) & set(cols_b_without_stain))
    
    # Merge the dataframes on the common columns
    merged_df = pd.merge(
        df_a, 
        df_b, 
        on=merge_on,
        how='outer'  # Use 'inner' if you only want rows present in both
    )
    
    # Save the merged dataframe
    merged_df.to_csv(args.output, index=False)
    print(f"Successfully merged {args.stainA_csv} and {args.stainB_csv} into {args.output}")
    print(f"Merged on columns: {merge_on}")


if __name__ == '__main__':
    main()