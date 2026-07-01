#!/usr/bin/env python3

import argparse
import numpy as np
import tifffile
import os

from imaris_ims_file_reader.ims import ims


def parse_args():
    parser = argparse.ArgumentParser(
        description="Subset channels from an IMS or TIFF image and save as ImageJ-compatible TIFF."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input file path (.ims or .tif/.tiff)"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output .tif file path"
    )

    parser.add_argument(
        "-c", "--channels",
        required=True,
        type=int,
        nargs="+",
        help="Zero-indexed channel indices to keep (e.g., -c 0 2 3)"
    )

    return parser.parse_args()


def find_channel_axis(shape):
    """
    Heuristic: channel axis is the smallest dimension.
    """
    return int(np.argmin(shape))


def subset_channels(data, channel_indices):
    """
    Subset the data along the detected channel axis.
    """
    channel_axis = find_channel_axis(data.shape)

    # Move channel axis to front
    data = np.moveaxis(data, channel_axis, 0)

    # Subset channels
    data = data[channel_indices]

    # Move axis back
    data = np.moveaxis(data, 0, channel_axis)

    return data


def load_image(path):
    """
    Load either IMS or TIFF automatically based on extension.
    """
    ext = os.path.splitext(path)[1].lower()

    if ext == ".ims":
        img = ims(path, squeeze_output=True)[0]
        data = np.array(img)
        return data

    elif ext in [".tif", ".tiff"]:
        return tifffile.imread(path)

    else:
        raise ValueError(f"Unsupported file type: {ext}")


def main():
    args = parse_args()

    # Load input image
    data = load_image(args.input)
    print("Image loaded shape: ", data.shape)
    # Subset channels
    subset = subset_channels(data, args.channels)

    # Save as ImageJ-compatible TIFF
    tifffile.imwrite(
        args.output,
        subset,
        imagej=True
    )

    print(f"Saved subset image to {args.output}")


if __name__ == "__main__":
    main()