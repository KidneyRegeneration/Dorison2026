import argparse
import tifffile
import numpy as np
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Find podocytes by combining DAPI and NPHS1 masks.')
    
    parser.add_argument('-n', '--nphs1_mask', 
                       required=True,
                       help='input nphs1 mask filepath')
    
    parser.add_argument('-d', '--dapi_mask',
                       required=True, 
                       help='input dapi filepath')
    
    parser.add_argument('-o', '--output',
                       required=True, 
                       help='output podocyte mask filepath')
    
    
    
    return parser.parse_args()


class MaskCombiner:
    def __init__(self, nphs1, dapi, output):
        self.nphs1 = self.read_mask_in(nphs1, binarise=True)
        self.dapi = self.read_mask_in(dapi)
        self.output = Path(output)
        
        
    def read_mask_in(self, mask_filepath, binarise=False):
        mask = tifffile.imread(mask_filepath)
        
        if binarise:
            mask = np.where(mask !=0, True, False)
            
        return mask

    def combine_masks(self):
        if self.nphs1.shape != self.dapi.shape:
            print(f"ERROR: MASK SHAPES DO NOT MATCH, NPHS1 SHAPE = {self.nphs1.shape}, DAPI SHAPE = {self.dapi.shape}")
            exit()
        else:
            print("Number of cells before masking: ", len(np.unique(self.dapi)))
            self.podocytes = np.where(self.nphs1, self.dapi, 0)
            print("Number of podocytes after masking: ", len(np.unique(self.podocytes)))
            
    def filter_wholly_contained_objects(self):
        unique_labels = np.unique(self.dapi)
        unique_labels = unique_labels[unique_labels != 0]
        
        # Create output mask
        self.podocytes = np.zeros_like(self.dapi)
        
        for label in unique_labels:
            # Get all pixels belonging to this object
            object_pixels = (self.dapi == label)
            
            # Check if ALL pixels of this object are within the binary mask
            if np.all(self.nphs1[object_pixels]):
                # Object is wholly contained, keep it
                self.podocytes[object_pixels] = label
        

    def save_out(self):
        tifffile.imwrite(self.output, self.podocytes)


def main():
    args = parse_args()
    combiner = MaskCombiner(args.nphs1_mask, args.dapi_mask, args.output)
    combiner.filter_wholly_contained_objects()
    combiner.save_out()
    

if __name__ == "__main__":
    main()