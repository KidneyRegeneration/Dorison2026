import argparse
import numpy as np
from imaris_ims_file_reader.ims import ims
import h5py
import tifffile
import pandas as pd

 
def parse_args():
    parser = argparse.ArgumentParser(description='Get stain amounts in each podocyte')
    
    parser.add_argument('-i', '--input_image',
                       required=True,
                       help='Input ims filepath')

    
    parser.add_argument('-p', '--pod_mask', 
                       required=True,
                       help='input image filepath')
    
    parser.add_argument('-d', '--nuclei_mask', 
                       required=True,
                       help='input nuclei mask filepath')
    
    parser.add_argument('-o', '--output_img',
                       required=True, 
                       type=str,
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
                       default=1,
                       type=int,
                       help='Index of stain (e.g. KAT2B) in image (zero-based indexing)')
    
    
    parser.add_argument('--extra_input',
                       required=False,
                       help='Extra input for nextflow rerun')
    
    
  
    return parser.parse_args()


class stainQuantifier:
    def __init__(self, input_image, pod_mask, output_mask, output_h5, output_csv, nuclei_mask, stain_index=1, dapi_index=0, ):
        self.input_img = self.read_image_in(input_image)
        self.pod_mask = self.read_mask_in(pod_mask)
        self.nuclei_mask = self.read_mask_in(nuclei_mask)
        self.output_img = output_mask
        self.output_h5 = output_h5
        self.stain_index = stain_index
        self.csv_out = output_csv
        self.dapi_index = dapi_index

    def read_image_in(self, filepath):
        image = ims(filepath)
        print(f"Image read in with shape = {image.shape}")
        return image

    def subset_stain(self):
        self.stain_channel = self.input_img[0, self.stain_index, 0, :, :]
        print(f"Image subset to stain chanel with shape = {self.stain_channel.shape}")
        self.check_shape()
        
    def subset_dapi(self):
        self.dapi_channel = self.input_img[0, self.dapi_index, 0, :, :]
        print(f"Image subset to stain chanel with shape = {self.dapi_channel.shape}")
        
    def read_mask_in(self, filepath):
        mask = tifffile.imread(filepath)
        print(f"Mask read in with shape = {mask.shape}")
        return mask
    
    def binarise_mask(self):
        self.binary_pod = np.where(self.pod_mask > 0, True, False)

    def check_shape(self):
        if self.stain_channel.shape != self.pod_mask.shape:
            print("ERROR: STAIN CHANNEL AND MASK SHAPE DO NOT MATCH!")
            exit()    
        
    def get_stain_averages(self):
        self.stain_mean = self.stain_channel.mean()
        self.stain_size = np.sum(np.where(self.stain_channel > 0, 1, 0))
        self.stain_mean_in_pods = self.stain_channel[self.pod_mask > 0].mean()
        self.stain_mean_outside_pods = self.stain_channel[self.pod_mask == 0].mean()
        self.stain_mean_outside_nuclei = self.stain_channel[self.nuclei_mask == 0].mean()
        self.stain_median = np.median(self.stain_channel)
        self.dapi_max = self.dapi_channel.max()
        self.dapi_min = self.dapi_channel.min()
        self.stain_max = self.stain_channel.max()
        self.stain_min = self.stain_channel.min()
        self.dapi_mean = self.dapi_channel.mean()
        
        print(f"Stats: Mean = {self.dapi_mean}, Median = {self.stain_median}, Max = {self.dapi_max}, Min = {self.dapi_min}")
    
    def min_max_norm(self):
        self.stain_minmax = ( self.stain_channel - self.dapi_min ) / (self.dapi_max - self.dapi_min)
        print(f"Min Max Normalised. Range was {self.stain_channel.min()} → {self.stain_channel.max()}, now: {self.stain_minmax.min()} → {self.stain_minmax.max()} ")
    
    def mask_stain(self):
        self.masked_stain = np.where(self.binary_pod, self.stain_channel, 0)
        self.masked_minmax = np.where(self.binary_pod, self.stain_minmax, 0)
    
    def get_stain_per_pod(self):
        self.stain_values = []
        self.pod_size = []
        self.stain_values_minmax = []
        self.pod_indices = []

        for pod_index in np.unique(self.pod_mask):
            if pod_index == 0:
                continue
            pod_mask = np.where(self.pod_mask == pod_index, True, False)
            stain_per_pod = pod_mask * self.stain_channel
            stain_per_pod_norm = pod_mask * self.stain_minmax
            self.stain_values.append(stain_per_pod.sum())
            self.pod_size.append(pod_mask.sum())
            self.stain_values_minmax.append(stain_per_pod_norm.sum())
            self.pod_indices.append(pod_index)

        self.norm_stain = np.array(self.stain_values) / np.array(self.pod_size)
        self.min_max_div_size = np.array(self.stain_values_minmax) / np.array(self.pod_size)
   
    def make_df(self):
        self.df = pd.DataFrame({"podocyte_index": self.pod_indices,
                                "stain_sum_in_pods": self.stain_values,
                                "stain_mean_overall": self.stain_mean,
                                "number_positive_stain_pixels": self.stain_size,
                                "pod_size": self.pod_size,
                                "stain_sum_norm_size": self.norm_stain,
                                "stain_sum_norm_minmax": self.stain_values_minmax,
                                "stain_sum_norm_min_max_norm_size": self.min_max_div_size,
                                "dapi_min": self.dapi_min,
                                "dapi_max": self.dapi_max,
                                "dapi_mean": self.dapi_mean,
                                "stain_mean_in_pods": self.stain_mean_in_pods,
                                "stain_mean_out_pods": self.stain_mean_outside_pods,
                                "stain_mean_out_nuclei": self.stain_mean_outside_nuclei,
        }
                               )
    
        self.df.to_csv(self.csv_out, index=False)
    
    def save_out(self):
        tifffile.imwrite(self.output_img, self.masked_stain)

        with h5py.File(self.output_h5, 'w') as f:
            # Save scalar variables
            f.create_dataset('img_mean', data=self.stain_mean)
            f.create_dataset('stain_mean_in_pods', data=self.stain_mean_in_pods)
            f.create_dataset('stain_mean_out_pods', data=self.stain_mean_outside_pods)
            f.create_dataset('stain_mean_out_nuclei', data=self.stain_mean_outside_nuclei)
            f.create_dataset('img_median', data=self.stain_median)
            f.create_dataset('img_max', data=self.stain_max)
            f.create_dataset('img_min', data=self.stain_min)
            f.create_dataset('dapi_max', data=self.dapi_max)
            f.create_dataset('dapi_min', data=self.dapi_min)
            f.create_dataset('dapi_mean', data=self.dapi_mean)
            
            # Convert to numpy arrays for storage
            f.create_dataset(f'stain_per_pod', data=self.stain_values)
            f.create_dataset(f'pod_sizes', data=self.pod_size)
            

    def process(self):
        self.binarise_mask()
        self.subset_stain()
        self.subset_dapi()
        self.get_stain_averages()
        self.min_max_norm()
        self.mask_stain()
        self.get_stain_per_pod()
        self.save_out()
        self.make_df()
        

def main():
    args = parse_args()
    stain_quantify = stainQuantifier(args.input_image, args.pod_mask, args.output_img, args.output_h5, args.output_csv, args.nuclei_mask, args.stain_index)
    stain_quantify.process()

if __name__ == "__main__":
    main()