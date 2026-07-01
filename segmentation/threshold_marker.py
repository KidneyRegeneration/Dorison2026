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
from skimage.filters import (
    threshold_otsu, threshold_triangle, threshold_li,
    threshold_yen, threshold_isodata, threshold_mean,
    threshold_minimum, threshold_niblack, threshold_sauvola
)
from skimage.restoration import rolling_ball
import argparse
from skimage.morphology import white_tophat, black_tophat, disk, ball

THRESHOLD_METHODS = {
    'otsu': threshold_otsu,
    'triangle': threshold_triangle,
    'li': threshold_li,
    'yen': threshold_yen,
    'isodata': threshold_isodata,
    'mean': threshold_mean,
    'minimum': threshold_minimum,
    'niblack': threshold_niblack,
    'sauvola': threshold_sauvola,
}

class MarkerFinder:
    def __init__(self, input_image, channel_index, output, threshold_method='otsu',
                 global_mask=None, mask_apply_stage='before',
                 rolling_ball_radius=None, smoothing_method="rolling_ball"):
        self.img = self.read_image_in(input_image)
        self.channel = self.img[0, channel_index, 0, :, :]
        print(f"Subsetting to channel index {channel_index}, shape is: {self.channel.shape}")
        self.output_image = Path(output)
        self.threshold_method = threshold_method
        self.mask_apply_stage = mask_apply_stage
        self.rolling_ball_radius = rolling_ball_radius
        self.smoothing_method = smoothing_method
        # Load global mask if provided
        self.global_mask = None
        if global_mask is not None:
            raw_mask = tifffile.imread(global_mask)
            raw_mask = np.squeeze(raw_mask)
            if raw_mask.shape != self.channel.shape:
                raise ValueError(
                    f"Global mask shape {raw_mask.shape} does not match "
                    f"channel shape {self.channel.shape}"
                )
            self.global_mask = (raw_mask > 0)
            print(f"Global mask loaded from {global_mask}, "
                  f"will be applied {mask_apply_stage} thresholding")

    def read_image_in(self, filepath):
        image = ims(filepath)
        print(f"Image read in with shape = {image.shape}")
        return image

    def _apply_mask_to_image(self, image):
        """Zero out pixels outside the global mask. Works on any dtype."""
        masked = image.copy()
        masked[~self.global_mask] = 0
        return masked

    def _apply_mask_to_binary(self, binary):
        """Restrict a binary mask to the global mask region."""
        return (binary & self.global_mask).astype(np.uint8)

    
    def _subtract_background(self, image):
        print("Subtracting background with rolling ball filter..")
        str_el = disk(self.rolling_ball_radius) 

        return white_tophat(image, str_el)


    def preprocess(self, save=False):
        channel = self.channel.copy()

        # Apply global mask before preprocessing if requested
        if self.global_mask is not None and self.mask_apply_stage == 'before':
            print("Applying global mask before preprocessing")
            channel = self._apply_mask_to_image(channel)
        img = channel
        # img = skimage.exposure.equalize_adapthist(channel, clip_limit=0.05)
        # img = skimage.exposure.adjust_gamma(img, 0.3)
        # self.pre_channel = skimage.exposure.adjust_sigmoid(img)
        # img = skimage.filters.gaussian(img, 2, preserve_range=True)
        img = self.smooth(img)
        # Background subtraction after smoothing, before normalisation
        # if self.rolling_ball_radius is not None:
        #     img = self._subtract_background(img)
        # Protect against an all-zero image after background subtraction
        img_max = img.max()
        if img_max == 0:
            print("Warning: image is all zeros after background subtraction — "
                  "check that rolling_ball_radius is not larger than your signal features.")
            self.preprocessed = img.astype(np.uint8)
        else:
            self.preprocessed = ((img / img_max) * 255).astype(np.uint8)

        if save:
            tifffile.imwrite(
                self.output_image.parents[0] / (self.output_image.stem + ".preprocessed" + self.output_image.suffix),
                self.preprocessed
            )

    def smooth(self, image_in):
        if self.smoothing_method == "gauss":
            img = skimage.filters.gaussian(image_in, 2, preserve_range=True)
        elif self.smoothing_method == "rolling_ball":
            img = self._subtract_background(image_in)
        else:
            img = image_in
        return img

    def threshold(self, save=False):
        gray = self.preprocessed.copy()

        thresh_fn = THRESHOLD_METHODS.get(self.threshold_method)
        if thresh_fn is None:
            raise ValueError(
                f"Unknown threshold method '{self.threshold_method}'. "
                f"Choose from: {list(THRESHOLD_METHODS.keys())}"
            )

        thresh = thresh_fn(gray)
        binary = (gray > thresh).astype(np.uint8)

        # Apply global mask to binary output after thresholding if requested
        if self.global_mask is not None and self.mask_apply_stage == 'after':
            print("Applying global mask after thresholding")
            binary = self._apply_mask_to_binary(binary)

        self.binary = binary

        if save:
            tifffile.imwrite(
                self.output_image.parents[0] / (self.output_image.stem + ".binary" + self.output_image.suffix),
                self.binary
            )

    def save_outputs(self):
        tifffile.imwrite(self.output_image, self.binary, imagej=True)

    def find_marker(self, save=False):
        self.preprocess(save)
        self.threshold(save)


def parse_args():
    parser = argparse.ArgumentParser(description='Threshold a single marker channel from an IMS image')

    parser.add_argument('-i', '--input',
                        required=True,
                        help='Input image filepath')

    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output mask filepath')

    parser.add_argument('-c', '--channel_index',
                        required=False,
                        default=0,
                        type=int,
                        help='Index of target channel in image (zero-based indexing)')

    parser.add_argument('-t', '--threshold_method',
                        required=False,
                        default='otsu',
                        choices=list(THRESHOLD_METHODS.keys()),
                        help='Thresholding method to use (default: otsu)')
    
    parser.add_argument('-s', '--smoothing_method',
                        required=False,
                        default='rolling-ball',
                        choices=["rolling_ball", "gauss"],
                        help='Thresholding method to use (default: otsu)')

    parser.add_argument('-g', '--global_mask',
                        required=False,
                        default=None,
                        help='Optional global mask TIFF to restrict analysis region')

    parser.add_argument('--mask_apply_stage',
                        required=False,
                        default='before',
                        choices=['before', 'after'],
                        help=(
                            'When to apply the global mask: '
                            '"before" zeros out the raw channel before preprocessing/thresholding '
                            '(affects threshold calculation); '
                            '"after" applies the mask to the binary result only '
                            '(threshold is computed on the full image). '
                            'Default: before'
                        ))

    parser.add_argument('--rolling_ball_radius',
                        required=False,
                        default=50,
                        type=float,
                        help=(
                            'Radius (in pixels) for rolling ball background subtraction. '
                            'Applied after Gaussian smoothing, before thresholding. '
                            'Should be set to roughly the radius of the largest background '
                            'feature you want to remove (i.e. larger than your signal structures). '
                            'Omit to skip background subtraction.'
                        ))

    parser.add_argument('--save_intermediates',
                        action='store_true',
                        help='Save intermediate masks for debugging')

    return parser.parse_args()


def main():
    args = parse_args()
    finder = MarkerFinder(
        args.input,
        args.channel_index,
        args.output,
        threshold_method=args.threshold_method,
        global_mask=args.global_mask,
        mask_apply_stage=args.mask_apply_stage,
        rolling_ball_radius=args.rolling_ball_radius,
        smoothing_method=args.smoothing_method
    )
    finder.find_marker(save=args.save_intermediates)
    finder.save_outputs()


if __name__ == "__main__":
    main()