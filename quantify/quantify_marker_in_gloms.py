import argparse
import numpy as np
from imaris_ims_file_reader.ims import ims
import h5py
import tifffile
import pandas as pd

 
def parse_args():
    parser = argparse.ArgumentParser(description='Get marker quantification in glomeruli')
    
    parser.add_argument('-i', '--input_image',
                       required=True,
                       help='Input ims filepath')

    parser.add_argument('-g', '--glom_mask', 
                       required=True,
                       help='Glom mask filepath')
    
    parser.add_argument('-o', '--output_img',
                       required=True, 
                       type=str,
                       nargs='+',
                       help='output mask filepath')
    
    parser.add_argument('-f', '--output_h5',
                       required=True,
                       type=str,
                       help='output h5 filepath')
    
    parser.add_argument('-c', '--output_csv',
                       required=True,
                       type=str,
                       help='output csv filepath')
    
    parser.add_argument('-n', '--stain_index',
                        required=False,
                        default=[1],
                        nargs='+',
                        type=int,
                        help='Index/indices of stain(s) in image (zero-based indexing). One or more space-separated integers.')

    parser.add_argument('--stain_names',
                        required=False,
                        nargs='+',
                        default=["marker"],
                        type=str,
                        help='Names of stains corresponding to each stain index. Must match the number of indices provided.')
    
    parser.add_argument('-x', '--nuclei_index',
                       required=False,
                       default=0,
                       type=int,
                       help='Index of nuclei marker (e.g. DAPI) in image (zero-based indexing)')
    
    
    parser.add_argument('--extra_input',
                       required=False,
                       help='Extra input for nextflow rerun')
    
    
  
    return parser.parse_args()


class stainQuantifier:
    def __init__(self, input_image, glom_mask, output_mask, output_h5, output_csv, stain_index=[1], dapi_index=0,
                 marker_names=['marker']):
        self.input_img = self.read_image_in(input_image)
        self.glom_mask = self.read_mask_in(glom_mask)
        self.output_img = output_mask
        self.output_h5 = output_h5
        self.stain_index = stain_index
        self.csv_out = output_csv
        self.dapi_index = dapi_index
        self.stain_names = marker_names

    def read_image_in(self, filepath):
        image = ims(filepath)
        print(f"Image read in with shape = {image.shape}")
        return image

    def subset_stain(self):
        self.stain_channel = []
        for i in self.stain_index:
            subset_stain_channel = self.input_img[0, i, 0, :, :]
            self.stain_channel.append(subset_stain_channel)
            print(f"Image subset to stain chanel with shape = {subset_stain_channel.shape}")
            self.check_shape(subset_stain_channel)
        
    def subset_dapi(self):
        self.dapi_channel = self.input_img[0, self.dapi_index, 0, :, :]
        print(f"Image subset to stain chanel with shape = {self.dapi_channel.shape}")
        
    def read_mask_in(self, filepath):
        mask = tifffile.imread(filepath)
        print(f"Mask read in with shape = {mask.shape}")
        return mask
    
    def binarise_mask(self):
        self.binary_glom = np.where(self.glom_mask > 0, True, False)

    def check_shape(self, img_to_check):
        if img_to_check.shape != self.glom_mask.shape:
            print("ERROR: STAIN CHANNEL AND MASK SHAPE DO NOT MATCH!")
            exit()    
            
    def get_stain_averages(self):
        self.dapi_max = self.dapi_channel.max()
        self.dapi_min = self.dapi_channel.min()
        self.dapi_mean = self.dapi_channel.mean()

        self.stain_mean = []
        self.stain_size = []
        self.stain_mean_in_gloms = []
        self.stain_mean_outside_gloms = []
        self.stain_median = []
        self.stain_max = []
        self.stain_min = []

        for channel in self.stain_channel:
            self.stain_mean.append(channel.mean())
            self.stain_size.append(np.sum(np.where(channel > 0, 1, 0)))
            self.stain_mean_in_gloms.append(channel[self.glom_mask > 0].mean())
            self.stain_mean_outside_gloms.append(channel[self.glom_mask == 0].mean())
            self.stain_median.append(np.median(channel))
            self.stain_max.append(channel[self.glom_mask > 0].max())
            self.stain_min.append(channel[self.glom_mask > 0].min())
            
            print(f"Stats: Mean = {self.dapi_mean}, Median = {self.stain_median}, Max = {self.dapi_max}, Min = {self.dapi_min}")
    
    def min_max_norm(self):
        self.stain_minmax = []
        for channel in self.stain_channel:
            minmaxed = (channel - self.dapi_min ) / (self.dapi_max - self.dapi_min)
            self.stain_minmax.append(minmaxed)
        print(f"Min Max Normalised. Range was {channel.min()} → {channel.max()}, now: {minmaxed.min()} → {minmaxed.max()} ")
    
    def mask_stain(self):
        self.masked_stain = []
        for channel in self.stain_channel:
            masked_stain = np.where(self.binary_glom, channel, 0)
            self.masked_stain.append(masked_stain)
    
    def get_stain_per_glom(self):
        self.stain_values = []
        self.glom_size = []
        self.stain_values_minmax = []
        self.norm_stain = []
        self.min_max_div_size = []
        
        for index, channel in enumerate(self.stain_channel):
            stain_channel_values = []
            stain_channel_minmax_values = []
            glom_sizes = []

            glom_mask = np.where(self.glom_mask > 0, True, False)
            stain_per_glom = glom_mask * channel
            stain_per_glom_norm = glom_mask * self.stain_minmax[index]
            stain_channel_values.append(stain_per_glom.sum())
            glom_sizes.append(glom_mask.sum())
            stain_channel_minmax_values.append(stain_per_glom_norm.sum())
            
            self.stain_values.append(stain_channel_values)
            self.glom_size.append(glom_sizes)
            self.stain_values_minmax.append(stain_channel_minmax_values)
            norm_stain = np.array(stain_channel_values) / np.array(glom_sizes)
            min_max_norm_size = np.array(stain_channel_minmax_values) / np.array(glom_sizes)       
            
            self.norm_stain.append(norm_stain)
            self.min_max_div_size.append(min_max_norm_size)
        
        print("length of stain values = ", len(self.stain_values))
        print("length of norm_stain values = ", len(self.norm_stain))
        print("length of min_max_div_size values = ", len(self.min_max_div_size))
        print("length of stain[0] values = ", len(self.stain_values[0]))
        print("length of glom[0] sizes = ", len(self.glom_size[0]))
   
    def make_df(self):

        base_df = pd.DataFrame({
            "glom_size": self.glom_size[0],
            "nuclei_channel_min": [self.dapi_min],
            "nuclei_channel_max": [self.dapi_max],
            "nuclei_channel_mean": [self.dapi_mean],
        })

        stain_dfs = []
        for i, name in enumerate(self.stain_names):
            stain_df = pd.DataFrame({
                f"{name}_total_marker_per_glom": self.stain_values[i],
                f"{name}_number_positive_stain_pixels": self.stain_size[i],
                f"{name}_total_marker_per_glom_divided_by_glom_size": self.norm_stain[i],
                f"{name}_total_marker_norm_by_minmax_of_nuclei_channel": self.stain_values_minmax[i],
                f"{name}_total_marker_norm_by_minmax_of_nuclei_channel_divided_by_glom_size": self.min_max_div_size[i],

                f"{name}_marker_channel_mean": [self.stain_mean[i]],
                f"{name}_mean_in_all_gloms": [self.stain_mean_in_gloms[i]],
                f"{name}_mean_outside_all_gloms": [self.stain_mean_outside_gloms[i]],
            })
            stain_dfs.append(stain_df)

        self.df = pd.concat([base_df] + stain_dfs, axis=1)
        self.df.to_csv(self.csv_out, index=False)
        
    def save_out(self):
        for index, i in enumerate(self.stain_channel):
            tifffile.imwrite(self.output_img[index], self.masked_stain[index])

        with h5py.File(self.output_h5, 'w') as f:
            # Save DAPI/nuclei scalars (unchanged)
            f.create_dataset('nuclei_channel_max', data=self.dapi_max)
            f.create_dataset('nuclei_channel_min', data=self.dapi_min)
            f.create_dataset('nuclei_channel_mean', data=self.dapi_mean)
            f.create_dataset('glom_sizes', data=self.glom_size)

            # Save per-stain datasets using stain name as key
            for i, name in enumerate(self.stain_names):
                f.create_dataset(f'{name}_img_mean', data=self.stain_mean[i])
                f.create_dataset(f'{name}_img_median', data=self.stain_median[i])
                f.create_dataset(f'{name}_marker_channel_max_in_gloms', data=self.stain_max[i])
                f.create_dataset(f'{name}_marker_channel_min_in_gloms', data=self.stain_min[i])
                f.create_dataset(f'{name}_mean_in_gloms', data=self.stain_mean_in_gloms[i])
                f.create_dataset(f'{name}_mean_out_gloms', data=self.stain_mean_outside_gloms[i])
                f.create_dataset(f'{name}_per_glom', data=self.stain_values[i])
            

    def process(self):
        self.binarise_mask()
        self.subset_stain()
        self.subset_dapi()
        self.get_stain_averages()
        self.min_max_norm()
        self.mask_stain()
        self.get_stain_per_glom()
        self.save_out()
        self.make_df()
        

def main():
    args = parse_args()
    if args.stain_names is not None and len(args.stain_names) != len(args.stain_index):
        print(f'Number of stain_names ({len(args.stain_names)}) must match '
                    f'number of stain_index values ({len(args.stain_index)})')
        exit()
    stain_quantify = stainQuantifier(args.input_image, args.glom_mask, args.output_img, args.output_h5, args.output_csv,
                                     args.stain_index, args.nuclei_index, args.stain_names)
    stain_quantify.process()

if __name__ == "__main__":
    main()