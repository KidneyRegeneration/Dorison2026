import argparse
import cv2
import numpy as np
from PIL import Image
from imaris_ims_file_reader.ims import ims
import tifffile
import skimage


def parse_args():
    parser = argparse.ArgumentParser(description='creates gif with image/mask overlay sequence')
    
    parser.add_argument('-i', '--image', 
                       required=True,
                       help='input image filepath')
    
    parser.add_argument('-m', '--mask',
                       required=True, 
                       help='input mask filepath')
    
    parser.add_argument('-o', '--output',
                       required=True,
                       help='output gif filepath')
    
    parser.add_argument('-n', '--nphs1_index',
                       required=False,
                       type=int,
                       default=2,
                       help='output gif filepath')
    
    return parser.parse_args()

def read_image_in(filepath, nphs1_index=2):
    image = ims(filepath)
    print(f"Image read in with shape = {image.shape}")
    nphs1 = image[0, nphs1_index, 0, :,  :]
    print(f"Subsetting to NPSH1 at index {nphs1_index}, shape is: {nphs1.shape} ")
    nphs1 = ((nphs1 / nphs1.max()) * 255).astype(np.uint8)
    return nphs1

def preprocess(img):
    img = skimage.exposure.equalize_adapthist(img, clip_limit=0.05)
    # img = skimage.exposure.adjust_gamma(img, 0.3)
    # img = skimage.exposure.adjust_sigmoid(img)
    # img = skimage.filters.gaussian(img, 5, preserve_range=True)
    img = ((img / img.max()) * 255 ).astype(np.uint8)
    
    return img

def read_mask(filepath):
    """placeholder for mask reading"""
    mask = tifffile.imread(filepath)
    mask = np.where(mask > 0, 1, 0) # ensure mask is binary
    return mask

def create_overlay(image, mask, alpha=0.5):
    """creates 50% opacity overlay of mask on image"""
    overlay = np.zeros_like(image)
    overlay[mask > 0] = [0, 255, 0] 
    
    return cv2.addWeighted(image, 1.0, overlay, alpha, 0)

def generate_gif(image_path, mask_path, output_path, nphs1_index):
    """generates 3-frame gif: image -> overlay -> mask"""
    
    # read inputs
    img = read_image_in(image_path, nphs1_index=nphs1_index)
    img = preprocess(img)
    mask = read_mask(mask_path)
    
    # resize mask to match image if needed
    if mask.shape[:2] != img.shape[:2]:
        mask = cv2.resize(mask, (img.shape[1], img.shape[0]), interpolation='INTER_NEAREST')
    
    # create frames
    frame1 = cv2.cvtColor(img, cv2.COLOR_GRAY2RGB)
    frame2 = create_overlay(frame1, mask, 0.2)
    frame3 = (np.repeat(np.expand_dims(mask, -1), 3, -1)) * 255
    # frames = [np.swapaxes(np.swapaxes(i, 0, -1), 1, 2) for i in [frame1, frame2, frame3]]
    frames = [frame1, frame2 ,frame3]

    # convert to PIL images
    pil_frames = [Image.fromarray(f.astype(np.uint8)) for f in frames]
    
    # save as gif
    pil_frames[0].save(
        output_path,
        save_all=True,
        append_images=pil_frames[1:],
        duration=800,  # ms per frame
        loop=0
    )


def main():
    args = parse_args()
    generate_gif(args.image, args.mask, args.output, args.nphs1_index)
    print(f"gif saved to {args.output}")

if __name__ == "__main__":
    main()