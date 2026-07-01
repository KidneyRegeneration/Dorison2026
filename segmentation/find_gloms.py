import pandas as pd
import os
import skimage.morphology
import numpy as np
import matplotlib.pyplot as plt
import h5py
import tifffile
from imaris_ims_file_reader.ims import ims
import cv2
from pathlib import Path
import skimage
from skimage.filters import threshold_otsu, sobel
import argparse
import scipy
from skimage.morphology import dilation, disk


class GlomFinder:
    def __init__(self, input_image, nphs1_index, output, nephrin_output=None, gauss_sigma=5, grow_region=False):
        self.img = self.read_image_in(input_image)
        self.nphs1 = self.img[0, nphs1_index, 0, :,  :]
        print(f"Subsetting to NPSH1 at index {nphs1_index}, shape is: {self.nphs1.shape} ")
        self.output_image = Path(output)
        self.nephrin_output = nephrin_output
        self.sigma = gauss_sigma
        self.grow_region = grow_region

    def read_image_in(self, filepath):
        image = ims(filepath)
        print(f"Image read in with shape = {image.shape}")
        return image
    
    def preprocess(self, exposure_fix=True, save=False):
        if exposure_fix:
            print("Exposure fix is TRUE.")
            img = skimage.exposure.equalize_adapthist(self.nphs1, clip_limit=0.05)
            img = skimage.exposure.adjust_gamma(img, 0.3)
            self.pre_nephrin = skimage.exposure.adjust_sigmoid(img)
        else: 
            print("Exposure fix is FALSE.")
            self.pre_nephrin = self.nphs1
        img = skimage.filters.gaussian(self.pre_nephrin,  self.sigma , preserve_range=True)
        self.preprocessed = ((img / img.max()) * 255 ).astype(np.uint8)
        
        if save:
            tifffile.imwrite(self.output_image.parents[0] / (self.output_image.stem + ".preprocessed" + self.output_image.suffix), self.preprocessed) 
    
    def threshold(self, save):
        # grayscale conversion
        gray = self.preprocessed.copy()
        
        # otsu threshold to create binary mask
        thresh = threshold_otsu(gray)
        binary = gray > thresh
        self.binary = binary.astype(np.uint8)
        self.binary = skimage.morphology.binary_closing(binary)
        if self.grow_region:
            self.binary = self.grow_thresholded_region()
        else:
            self.binary = self.binary.astype(np.uint8)
        if save:
            tifffile.imwrite(self.output_image.parents[0] / (self.output_image.stem + ".binary" + self.output_image.suffix), self.binary) 
    
    
    def grow_thresholded_region(self): 
        grown = self.binary.copy()
        for _ in range(20):  # number of dilations
            dilated = dilation(grown, disk(3))
        grown = np.where(dilated, 1, 0).astype(np.uint8)
        return grown

    def threshold_for_nephrin(self):
        # otsu threshold to create binary mask
        thresh_neph = threshold_otsu(self.pre_nephrin)
        binary = self.pre_nephrin > thresh_neph
        self.binary_neph = binary.astype(np.uint8)
 
    def get_filled_contours(self):
        # find contours of filled mask
        self.contours, _ = cv2.findContours(self.binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    def mask_inside_contours(self):
        # create blank mask same size as image
        mask = np.zeros(self.preprocessed.shape[:2], dtype=np.uint8)
        
        # fill inside all contours
        cv2.fillPoly(mask, self.contours, 255)
        
        mask = skimage.morphology.diameter_opening(mask, 100)
        mask = skimage.morphology.dilation(mask, footprint=[(np.ones((9, 1)), 1), (np.ones((1, 9)), 1)], mode='nearest')
        mask = scipy.ndimage.binary_fill_holes(mask)

        self.masked_gloms = mask.astype(np.uint8)
        
    def save_outputs(self):
        tifffile.imwrite(self.output_image, self.masked_gloms, imagej=True)

    def _remove_edge_effect(self, mask, edge_percent=0.05):
        """
        Zero out edge pixels from all sides of a binary mask.
        
        Parameters:
        -----------
        mask : numpy.ndarray
            Binary mask array
        edge_percent : float
            Percentage (as decimal) of height/width to zero out from each side
            
        Returns:
        --------
        numpy.ndarray
            Modified mask with edges zeroed out
        """
        # Create a copy to avoid modifying the original
        result = mask.copy()
        
        height, width = mask.shape[:2]
        
        # Calculate number of pixels to remove from each side
        top_bottom_margin = int(height * edge_percent)
        left_right_margin = int(width * edge_percent)
        
        # Zero out top and bottom
        result[:top_bottom_margin, :] = 0
        result[-top_bottom_margin:, :] = 0
        
        # Zero out left and right
        result[:, :left_right_margin] = 0
        result[:, -left_right_margin:] = 0
        
        return result


    def process_nephrin_mask(self):
        self.threshold_for_nephrin()
        edgeless = self._remove_edge_effect(self.binary_neph)
        boolean_mask = np.where(edgeless > 0, True, False)
        # dilated_mask = skimage.morphology.dilation(self.binary, footprint=[(np.ones((9, 1)), 1), (np.ones((1, 9)), 1)], mode='nearest')
        # label_img = skimage.morphology.label(dilated_mask)
        nephrin_mask = skimage.morphology.remove_small_objects(boolean_mask, 400)
        # nephrin_mask = skimage.morphology.remove_objects_by_distance(label_img, boolean_mask.shape[1] // 4)
        # nephrin_mask = boolean_mask * nephrin_mask
        self.nephrin_mask = nephrin_mask.astype(np.uint8)

        
    def save_nephrin_mask(self):
        tifffile.imwrite(self.nephrin_output, self.nephrin_mask, imagej=True)

    
    def find_gloms(self, adjust_exposure, save=False):
        self.preprocess(exposure_fix=adjust_exposure, save=save)
        self.threshold(save)
        if self.nephrin_output is not None:
            self.process_nephrin_mask()
            self.save_nephrin_mask()
        self.get_filled_contours()
        self.mask_inside_contours()
        
        
def parse_args():
    parser = argparse.ArgumentParser(description='find glomeruli by thresholding NPHS1 channel')
    
    parser.add_argument('-i', '--input', 
                       required=True,
                       help='input image filepath')
    
    parser.add_argument('-o', '--output',
                       required=True, 
                       help='output mask filepath')
    
    parser.add_argument('-n', '--nphs1_index',
                       required=False,
                       default=2,
                       type=int,
                       help='Index of NPHS1 in image (zero-based indexing)')
    parser.add_argument('--nephrin_output',
                        required=False,
                       default=None,
                       help='Output for nephrin only mask')
    
    parser.add_argument('--save_intermediates',
                       action="store_true",
                       help='Save intermediate masks for debug')
    
    parser.add_argument('--nextflow_extra',
                       action="store_true",
                       help='Save intermediate masks for debug')
    
    parser.add_argument('-s', '--smooth_sigma',
                    required=False,
                    default=5,
                    type=int,
                    help='Gaussian blur sigma')
    
    parser.add_argument('-e', '--exposure_adjust',
                    required=False,
                    action="store_false",
                    help='TURN OFF adjust exposure before masking.')
    
    parser.add_argument('-g', '--grow_region',
                    required=False,
                    action="store_true",
                    help='Grow region after masking.')
    
    return parser.parse_args()

def main():
    args = parse_args()
    glom_finder = GlomFinder(args.input, args.nphs1_index, args.output, args.nephrin_output, args.smooth_sigma, args.grow_region)
    glom_finder.find_gloms(adjust_exposure=args.exposure_adjust, save=args.save_intermediates )
    glom_finder.save_outputs()
    


if __name__ == "__main__":
    main()