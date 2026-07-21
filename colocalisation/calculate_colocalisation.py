#!/usr/bin/env python3
"""
colocalize.py — Fluorescent protein colocalisation analysis.

Computes Modified Manders Coefficients (M1: A-in-B, M2: B-in-A) and Pearson's
Correlation Coefficient (PCC) for two fluorescent markers, and produces a scatter
plot and a representative snapshot of the colocalised region.

Modified Manders:
    M1 (A in B) = sum(A pixels where B_mask > 0) / sum(all A pixels)
    M2 (B in A) = sum(B pixels where A_mask > 0) / sum(all B pixels)

Supported formats
-----------------
  * TIFF / PNG / standard image formats (via Pillow / tifffile)
  * Imaris IMS files (via imaris_ims_file_reader)
    IMS volumes have 5 axes: (t, c, z, y, x). Use --timepoint and
    --z-slice to select the 2-D plane to analyse (defaults: t=0, z=0).

Input modes
-----------
  Two separate single-channel images with separate masks:
      colocalize.py -i chA.tif chB.tif --mask-a maskA.tif --mask-b maskB.tif -o results.csv

  Single multi-channel image with explicit channel indices:
      colocalize.py -i stack.tif --channels 0 2 --mask-a maskA.tif --mask-b maskB.tif -o results.csv

  With optional global mask to restrict analysis region:
      colocalize.py -i chA.tif chB.tif --mask-a mA.tif --mask-b mB.tif --mask mG.tif -o results.csv

  IMS file, channels 0 and 1, z-slice 5:
      colocalize.py -i dataset.ims --channels 0 1 --mask-a mA.tif --mask-b mB.tif --z-slice 5 -o results.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from skimage.morphology import footprint_rectangle
from skimage.filters import (
    threshold_otsu, threshold_triangle, threshold_li,
    threshold_yen, threshold_isodata, threshold_mean,
    threshold_minimum, threshold_niblack, threshold_sauvola
)
import skimage
try:
    import tifffile
    _TIFFFILE_AVAILABLE = True
except ImportError:
    _TIFFFILE_AVAILABLE = False

try:
    from imaris_ims_file_reader.ims import ims as ImsReader
    _IMS_AVAILABLE = True
except ImportError:
    _ImsReader = None
    _IMS_AVAILABLE = False


# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------

_IMS_SUFFIXES = {".ims"}
_TIFF_SUFFIXES = {".tif", ".tiff"}


def _is_ims(path: Path) -> bool:
    return path.suffix.lower() in _IMS_SUFFIXES


def _is_tiff(path: Path) -> bool:
    return path.suffix.lower() in _TIFF_SUFFIXES


# ---------------------------------------------------------------------------
# I/O helpers — IMS
# ---------------------------------------------------------------------------

def _open_ims(path: Path) -> "ImsReader":
    if not _IMS_AVAILABLE:
        raise RuntimeError(
            f"Cannot open '{path}': imaris_ims_file_reader is not installed. "
            "Install it with: pip install imaris-ims-file-reader"
        )
    return ImsReader(str(path))


def _validate_ims_axes(
    img: "ImsReader",
    channel: int,
    timepoint: int,
    z_slice: int,
    path: Path,
) -> None:
    shape = img.shape
    try:
        n_t, n_c, n_z = shape[0], shape[1], shape[2]
    except (TypeError, IndexError):
        return

    if timepoint < 0 or timepoint >= n_t:
        raise ValueError(f"Timepoint {timepoint} out of range (0–{n_t - 1}) for '{path}'.")
    if channel < 0 or channel >= n_c:
        raise ValueError(f"Channel {channel} out of range (0–{n_c - 1}) for '{path}'.")
    if z_slice < 0 or z_slice >= n_z:
        raise ValueError(f"Z-slice {z_slice} out of range (0–{n_z - 1}) for '{path}'.")


def load_ims_plane(path: Path, channel: int, timepoint: int, z_slice: int) -> np.ndarray:
    img = _open_ims(path)
    _validate_ims_axes(img, channel, timepoint, z_slice, path)
    return np.array(img[timepoint, channel, z_slice, :, :], dtype=np.float32)


def load_ims_single_channel_plane(path: Path, timepoint: int, z_slice: int) -> np.ndarray:
    img = _open_ims(path)
    try:
        n_c = img.shape[1]
    except (TypeError, IndexError):
        n_c = None

    if n_c is not None and n_c != 1:
        raise ValueError(
            f"'{path}' has {n_c} channels. Use --channels to select which "
            "two channels to compare, or provide two single-channel IMS files."
        )
    channel = 0
    _validate_ims_axes(img, channel, timepoint, z_slice, path)
    return np.array(img[timepoint, channel, z_slice, :, :], dtype=np.float32)


# ---------------------------------------------------------------------------
# I/O helpers — TIFF / standard images
# ---------------------------------------------------------------------------

def _load_tiff_array(path: Path) -> np.ndarray:
    if _TIFFFILE_AVAILABLE:
        arr = tifffile.imread(str(path)).astype(np.float32)
    else:
        arr = np.array(Image.open(path), dtype=np.float32)
    return arr


def load_image_2d(path: Path) -> np.ndarray:
    if _is_tiff(path):
        arr = _load_tiff_array(path)
    else:
        arr = np.array(Image.open(path), dtype=np.float32)

    arr = np.squeeze(arr)

    if arr.ndim != 2:
        raise ValueError(
            f"'{path}' has shape {arr.shape} after squeezing. "
            "It appears to be multi-channel — use --channels to select channels."
        )
    return arr


def load_channel(path: Path, channel: int) -> np.ndarray:
    if _is_tiff(path):
        arr = _load_tiff_array(path)
    else:
        arr = np.array(Image.open(path), dtype=np.float32)

    arr = np.squeeze(arr)

    if arr.ndim == 2:
        if channel != 0:
            raise ValueError(
                f"'{path}' is a single-plane (2-D) image; channel index must be 0, got {channel}."
            )
        return arr

    if arr.ndim != 3:
        raise ValueError(f"Unexpected array shape {arr.shape} for '{path}'. Expected 2-D or 3-D.")

    n_channels = arr.shape[0]
    if channel < 0 or channel >= n_channels:
        raise ValueError(
            f"Channel index {channel} out of range for image with {n_channels} channel(s) at '{path}'."
        )
    return arr[channel]


def load_mask(path: Path, shape: tuple[int, int]) -> np.ndarray:
    """
    Load a binary mask and verify it matches expected spatial (H, W) shape.
    Returns a boolean array.
    """
    if _is_tiff(path):
        arr = _load_tiff_array(path).astype(bool)
    else:
        arr = np.array(Image.open(path), dtype=bool)

    arr = np.squeeze(arr)

    if arr.ndim == 3:
        if arr.shape[0] <= arr.shape[-1]:
            arr = arr.any(axis=0)   # (C, H, W) -> (H, W)
        else:
            arr = arr.any(axis=-1)  # (H, W, C) -> (H, W)

    if arr.shape != shape:
        raise ValueError(f"Mask shape {arr.shape} does not match image shape {shape}.")
    return arr


def smooth(image, smoothing_method):
    if smoothing_method == "gauss":
        print("Smoothing with gaussian blur")
        img = skimage.filters.gaussian(image, 2, preserve_range=True)
    elif smoothing_method == "rolling_ball":
        print("Smoothing with rolling ball")
        img = skimage.filters.rank.mean(image.astype(np.uint16), footprint=footprint_rectangle((3,3)))
        bg = skimage.restoration.rolling_ball(img, radius=50)
        img = img - bg
    else:
        print("SMOOTHING METHOD IS NOT CURRENTLY SUPPORTED. RETURNING ORIGINAL IMAGE.")
        img = image
    return img

def threshold(image, method):
    threshold_methods = {
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
    
    thresh_fn = threshold_methods.get(method)
    if thresh_fn is None:
        raise ValueError(
            f"Unknown threshold method '{method}'. "
            f"Choose from: {list(threshold_methods.keys())}"
        )
    thresh = thresh_fn(image)
    thresholded_image  = np.where(image > thresh, image, 0)
    return thresholded_image
    
def preprocess(image, smoothing_algo=None, thresholding_algo=None):
    if smoothing_algo is not None:
        print("Smoothing!")
        img = smooth(image, smoothing_algo)
    else:
        img = image
        
    if thresholding_algo is not None:
        print("Thresholding!")
        img = threshold(img, thresholding_algo)
    else:
        img = img
    
    return img

def modified_manders(
    signal: np.ndarray,
    coloc_mask: np.ndarray,
    global_mask: np.ndarray | None = None,
    smooth=None,
    thresh=None
) -> float:
    """
    Modified Manders coefficient for *signal* (A) restricted to *coloc_mask* (B).

        Without a global mask G:
            M(A in B) = sum(A * [B > 0]) / sum(A)

        With a global mask G:
            M(A in B) = sum(A * [B > 0] * [G > 0]) / sum(A * [G > 0])

        where [.] is 1 where the condition holds and 0 otherwise, and sums run
        over all pixels.

    Parameters
    ----------
    signal      : 2-D float array, the raw intensity channel being measured (A).
    coloc_mask  : boolean 2-D array, positive pixels of the *other* marker (B).
    global_mask : optional boolean 2-D array; when provided, both the
                  numerator and denominator are restricted to pixels within
                  this mask (i.e. coloc_mask is also intersected with it).

    Returns NaN when total signal within the global mask is zero.
    """
    print("Signal shape = ", signal.shape)
    if smooth is not None or thresh is not None:
        print("Preprocessing the signal with smoothing = ", smooth, "and thresholding = ", thresh)
        signal_out = preprocess(signal, smoothing_algo=smooth, thresholding_algo=thresh)
        signal = np.where(signal_out > 0, signal, 0) # test thresholding the signal but retain original values.

    if global_mask is not None:
        denom_pixels = np.where(global_mask > 0, signal, 0)
        numer_pixels = np.where(coloc_mask > 0, signal, 0)
        numer_pixels = np.where(global_mask > 0, numer_pixels, 0)
    else:
        denom_pixels = signal
        numer_pixels = np.where(coloc_mask > 0, signal, 0)
    
    total = denom_pixels.sum()
    if total == 0.0:
        return float("nan")
    numer_total = numer_pixels.sum()

    return float(numer_total / total)


def pearson_correlation_coefficient(
    a: np.ndarray,
    b: np.ndarray,
) -> float:
    """
    Pearson's Correlation Coefficient (PCC) between pixel intensities.
    Accepts either 2-D arrays (full image) or 1-D ravelled pixel vectors.
    Returns NaN when either channel has zero variance.
    """
    a_flat = a.ravel()
    b_flat = b.ravel()
    if a_flat.std() == 0.0 or b_flat.std() == 0.0:
        return float("nan")
    r, _ = pearsonr(a_flat, b_flat)
    return float(r)


def save_scatter_plot(
    a: np.ndarray,
    b: np.ndarray,
    name_a: str,
    name_b: str,
    output_path: Path,
    pcc: float,
    moc_a_in_b: float,
    moc_b_in_a: float,
    max_points: int = 50_000,
) -> None:
    rng = np.random.default_rng(seed=42)
    flat_a = a.ravel()
    flat_b = b.ravel()

    if flat_a.size > max_points:
        idx = rng.choice(flat_a.size, size=max_points, replace=False)
        flat_a = flat_a[idx]
        flat_b = flat_b[idx]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(flat_a, flat_b, s=1, alpha=0.3, color="steelblue", rasterized=True)
    ax.set_xlabel(f"{name_a} intensity", fontsize=12)
    ax.set_ylabel(f"{name_b} intensity", fontsize=12)
    ax.set_title(f"Colocalisation: {name_a} vs {name_b}", fontsize=13)
    annotation = (
        f"PCC           = {pcc:.4f}\n"
        f"MOC ({name_a} in {name_b}) = {moc_a_in_b:.4f}\n"
        f"MOC ({name_b} in {name_a}) = {moc_b_in_a:.4f}"
    )
    ax.text(
        0.05, 0.95, annotation,
        transform=ax.transAxes,
        va="top", ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7),
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"  Scatter plot saved → {output_path}")


# ---------------------------------------------------------------------------
# Output: representative snapshot
# ---------------------------------------------------------------------------

def _local_correlation(a: np.ndarray, b: np.ndarray, tile_size: int) -> tuple[int, int]:
    h, w = a.shape
    best_r = -np.inf
    best_row, best_col = 0, 0

    for row in range(0, h - tile_size + 1, tile_size):
        for col in range(0, w - tile_size + 1, tile_size):
            tile_a = a[row: row + tile_size, col: col + tile_size].ravel()
            tile_b = b[row: row + tile_size, col: col + tile_size].ravel()
            if tile_a.std() == 0 or tile_b.std() == 0:
                continue
            r, _ = pearsonr(tile_a, tile_b)
            if r > best_r:
                best_r = r
                best_row, best_col = row, col

    return best_row, best_col


def _normalise_to_uint8(arr: np.ndarray) -> np.ndarray:
    lo, hi = arr.min(), arr.max()
    if hi == lo:
        return np.zeros_like(arr, dtype=np.uint8)
    return ((arr - lo) / (hi - lo) * 255.0).astype(np.uint8)


def save_snapshot(
    a: np.ndarray,
    b: np.ndarray,
    name_a: str,
    name_b: str,
    output_path: Path,
    tile_size: int = 128,
) -> None:
    h, w = a.shape
    effective_tile = min(tile_size, h, w)
    row, col = _local_correlation(a, b, effective_tile)

    tile_a = a[row: row + effective_tile, col: col + effective_tile]
    tile_b = b[row: row + effective_tile, col: col + effective_tile]

    r_ch = _normalise_to_uint8(tile_a)
    g_ch = _normalise_to_uint8(tile_b)
    blue = np.zeros_like(r_ch)

    rgb = np.stack([r_ch, g_ch, blue], axis=-1)
    Image.fromarray(rgb, mode="RGB").save(output_path)
    print(
        f"  Snapshot saved → {output_path}  "
        f"(tile [{row}:{row+effective_tile}, {col}:{col+effective_tile}], "
        f"{name_a}=red, {name_b}=green)"
    )


def calculate_marker_in_gloms(marker_img, glom_mask):
    marker_subset_to_glom_mask = np.where(glom_mask > 0, marker_img, 0)
    total_marker_in_glom = marker_subset_to_glom_mask.sum()
    glom_size = np.count_nonzero(glom_mask)
    marker_in_glom_divided_by_glom_size = total_marker_in_glom / glom_size
    
    return marker_in_glom_divided_by_glom_size


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="colocalize.py",
        description=(
            "Compute modified Manders coefficients (MOC_AinB, MOC_BinA) and "
            "Pearson's Correlation Coefficient (PCC) for two fluorescent markers. "
            "Requires per-marker binary masks to define each marker's positive region."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # --- Input: raw images ---
    input_group = parser.add_argument_group("Input — raw images")
    input_group.add_argument(
        "-i", "--images",
        nargs="+",
        required=True,
        metavar="IMAGE",
        help=(
            "One or two image file paths (.tif/.tiff or .ims). "
            "Provide TWO single-channel images (one per marker), or ONE "
            "multi-channel image together with --channels."
        ),
    )
    input_group.add_argument(
        "-c", "--channels",
        nargs=2,
        type=int,
        metavar=("CH_A", "CH_B"),
        default=None,
        help=(
            "Zero-based channel indices to extract from a single multi-channel "
            "image. Required when --images receives exactly one file."
        ),
    )
    input_group.add_argument(
        "-z", "--z-slice",
        type=int,
        default=0,
        metavar="Z",
        help="Z-slice index (IMS / multi-Z TIFF). Default: 0.",
    )
    input_group.add_argument(
        "-t", "--timepoint",
        type=int,
        default=0,
        metavar="T",
        help="Time-point index (IMS only). Default: 0.",
    )

    # --- Input: masks ---
    mask_group = parser.add_argument_group("Input — masks")
    mask_group.add_argument(
        "--mask-a",
        required=True,
        metavar="MASK_A",
        help=(
            "Binary mask defining positive pixels for marker A. "
            "Used as the colocalisation region when computing MOC_BinA, "
            "and as the denominator region is always total signal within "
            "the optional global mask."
        ),
    )
    mask_group.add_argument(
        "--mask-b",
        required=True,
        metavar="MASK_B",
        help="Binary mask defining positive pixels for marker B.",
    )
    mask_group.add_argument(
        "-m", "--mask",
        metavar="GLOBAL_MASK",
        default=None,
        help=(
            "Optional global binary mask (same H × W as marker images). "
            "When provided, all analysis — numerator, denominator, and PCC — "
            "is restricted to pixels within this mask."
        ),
    )

    # --- Input: names ---
    input_group.add_argument(
        "-n", "--names",
        nargs=2,
        metavar=("NAME_A", "NAME_B"),
        default=None,
        help="Human-readable marker names. Defaults to 'markerA' and 'markerB'.",
    )
    
    input_group.add_argument(
        "--smooth_a",
        default=None,
        help="Smoothing algorithm to use for image A."
        
    )
    input_group.add_argument(
        "--smooth_b",
        default=None,
        help="Smoothing algorithm to use for image B."
        
    )
    
    input_group.add_argument(
        "--thresh_a",
        default=None,
        help="Thresholding algorithm to use for image A."
        
    )
    input_group.add_argument(
        "--thresh_b",
        default=None,
        help="Thresholding algorithm to use for image B."

    )
    
    # --- Output ---
    output_group = parser.add_argument_group("Output")
    output_group.add_argument(
        "-o", "--output-csv",
        required=True,
        metavar="CSV",
        help="Destination CSV file for colocalisation coefficients.",
    )
    output_group.add_argument(
        "-p", "--output-plot",
        metavar="PLOT",
        default=None,
        help="Scatter plot output path. Inferred from --output-csv when omitted.",
    )
    output_group.add_argument(
        "-s", "--output-snapshot",
        metavar="SNAPSHOT",
        default=None,
        help="Representative snapshot PNG path. Inferred from --output-csv when omitted.",
    )
    output_group.add_argument(
        "--tile-size",
        type=int,
        default=128,
        metavar="N",
        help="Side length (pixels) of the representative snapshot tile. Default: 128.",
    )
    

    return parser


def validate_args(args: argparse.Namespace) -> None:
    n_images = len(args.images)
    if n_images == 1 and args.channels is None:
        print("error: when providing a single image, --channels CH_A CH_B is required.", file=sys.stderr)
        sys.exit(1)
    if n_images == 2 and args.channels is not None:
        print("error: --channels is only valid when providing a single image.", file=sys.stderr)
        sys.exit(1)
    if n_images > 2:
        print(f"error: at most two images may be provided; got {n_images}.", file=sys.stderr)
        sys.exit(1)
    if args.z_slice < 0:
        print("error: --z-slice must be >= 0.", file=sys.stderr)
        sys.exit(1)
    if args.timepoint < 0:
        print("error: --timepoint must be >= 0.", file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args)

    # --- Resolve output paths ---
    csv_path = Path(args.output_csv)
    stem = csv_path.stem
    plot_path = Path(args.output_plot) if args.output_plot else csv_path.with_name(f"{stem}_scatter.png")
    snapshot_path = Path(args.output_snapshot) if args.output_snapshot else csv_path.with_name(f"{stem}_snapshot.png")

    name_a, name_b = args.names if args.names else ("markerA", "markerB")

    # --- Load raw intensity images ---
    print("Loading images …")
    t, z = args.timepoint, args.z_slice

    if len(args.images) == 2:
        path_a, path_b = Path(args.images[0]), Path(args.images[1])
        ch_a = (
            load_ims_single_channel_plane(path_a, timepoint=t, z_slice=z)
            if _is_ims(path_a) else load_image_2d(path_a)
        )
        ch_b = (
            load_ims_single_channel_plane(path_b, timepoint=t, z_slice=z)
            if _is_ims(path_b) else load_image_2d(path_b)
        )
    else:
        img_path = Path(args.images[0])
        ch_idx_a, ch_idx_b = args.channels
        if _is_ims(img_path):
            ch_a = load_ims_plane(img_path, channel=ch_idx_a, timepoint=t, z_slice=z)
            ch_b = load_ims_plane(img_path, channel=ch_idx_b, timepoint=t, z_slice=z)
        else:
            ch_a = load_channel(img_path, ch_idx_a)
            ch_b = load_channel(img_path, ch_idx_b)

    if ch_a.shape != ch_b.shape:
        print(f"error: marker image shapes do not match: {ch_a.shape} vs {ch_b.shape}.", file=sys.stderr)
        sys.exit(1)

    spatial_shape: tuple[int, int] = ch_a.shape

    # --- Load per-marker masks ---
    print("Loading marker masks …")
    mask_a = load_mask(Path(args.mask_a), spatial_shape)
    mask_b = load_mask(Path(args.mask_b), spatial_shape)

    # --- Load optional global mask ---
    global_mask: np.ndarray | None = None
    if args.mask:
        print("Loading global mask …")
        global_mask = load_mask(Path(args.mask), spatial_shape)

    # --- Compute modified Manders coefficients ---
    # MOC_AinB: fraction of A signal that overlaps with B-positive pixels
    # MOC_BinA: fraction of B signal that overlaps with A-positive pixels
    print("Computing colocalisation coefficients …")

    moc_a_in_b = modified_manders(ch_a, coloc_mask=mask_b, global_mask=global_mask, smooth=args.smooth_a, thresh=args.thresh_a)
    moc_b_in_a = modified_manders(ch_b, coloc_mask=mask_a, global_mask=global_mask, smooth=args.smooth_b, thresh=args.thresh_b)
    moc_a_in_a = modified_manders(ch_a, coloc_mask=mask_a, global_mask=global_mask, smooth=args.smooth_a, thresh=args.thresh_a)
    moc_b_in_b = modified_manders(ch_b, coloc_mask=mask_b, global_mask=global_mask, smooth=args.smooth_b, thresh=args.thresh_b)
    # --- PCC: computed over global mask region (or full image) ---
    if global_mask is not None:
        pcc_a = ch_a[global_mask]
        pcc_b = ch_b[global_mask]
    else:
        pcc_a = ch_a.ravel()
        pcc_b = ch_b.ravel()

    pcc = pearson_correlation_coefficient(pcc_a, pcc_b)

    mean_a = float(pcc_a.mean())
    mean_b = float(pcc_b.mean())

    markerA_in_gloms = calculate_marker_in_gloms(ch_a, global_mask)
    markerB_in_gloms = calculate_marker_in_gloms(ch_b, global_mask)


    print(f"  {name_a} mean intensity  : {mean_a:.4f}")
    print(f"  {name_b} mean intensity  : {mean_b:.4f}")
    print(f"  MOC ({name_a} in {name_b}) : {moc_a_in_b:.6f}")
    print(f"  MOC ({name_b} in {name_a}) : {moc_b_in_a:.6f}")
    print(f"  MOC ({name_a} in {name_a}) : {moc_a_in_a:.6f}")
    print(f"  MOC ({name_b} in {name_b}) : {moc_b_in_b:.6f}")
    print(f"  PCC                   : {pcc:.6f}")

    # --- Save CSV ---
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame([{
        f"MOC_{name_a}in{name_b}":      moc_a_in_b,
        f"MOC_{name_b}in{name_a}":      moc_b_in_a,
        f"MOC_{name_a}in{name_a}":      moc_a_in_a,
        f"MOC_{name_b}in{name_b}":      moc_b_in_b,
        "PCC":           pcc,
        f"{name_a}_mean":  mean_a,
        f"{name_b}_mean":  mean_b,
        f"{name_a}_in_gloms_normalised": markerA_in_gloms,
        f"{name_b}_in_gloms_normalised": markerB_in_gloms
    }])
    df.to_csv(csv_path, index=False)
    print(f"  Results saved → {csv_path}")

    # --- Scatter plot ---
    print("Generating scatter plot …")
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    save_scatter_plot(
        a=pcc_a, b=pcc_b,
        name_a=name_a, name_b=name_b,
        output_path=plot_path,
        pcc=pcc,
        moc_a_in_b=moc_a_in_b,
        moc_b_in_a=moc_b_in_a,
    )

    # --- Representative snapshot ---
    print("Generating representative snapshot …")
    snapshot_path.parent.mkdir(parents=True, exist_ok=True)
    snap_a = np.where(global_mask, ch_a, 0.0) if global_mask is not None else ch_a
    snap_b = np.where(global_mask, ch_b, 0.0) if global_mask is not None else ch_b
    save_snapshot(
        a=snap_a, b=snap_b,
        name_a=name_a, name_b=name_b,
        output_path=snapshot_path,
        tile_size=args.tile_size,
    )

    print("Done.")


if __name__ == "__main__":
    main()