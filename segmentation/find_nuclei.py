import numpy as np
import matplotlib.pyplot as plt
import tifffile
from imaris_ims_file_reader.ims import ims
import argparse
from cellpose import models, plot
import torch


class NucleiSegmentor:
    def __init__(self, image, dapi_index=0):
        print("Reading in...")
        self.image = ims(image)
        print("Image shape = ", self.image.shape)
        self.dapi_channel = self.image[0, dapi_index, :]
        print("DAPI channel shape = ", self.dapi_channel.shape)
    
    def find_nuclei_cellpose(self, diameter):
        gpu = torch.cuda.is_available()

        self.model = models.CellposeModel(gpu=gpu)
        flow_threshold = 0.4
        cellprob_threshold = 0.0
        tile_norm_blocksize = 0
        self.masks, self.flows, styles = self.model.eval(self.dapi_channel,
                                                         diameter=diameter,
                                                         batch_size=32,
                                                         flow_threshold=flow_threshold,
                                                         cellprob_threshold=cellprob_threshold,
                                                         normalize={"tile_norm_blocksize": tile_norm_blocksize})
        

    def save_outputs(self, output, binary_output=None):
        if not output.endswith(".tif"):
            print("WARNING: You're not saving the output mask as a tif, you probably should.")
        tifffile.imwrite(output, self.masks)
        print("Saved uniquely labelled mask to: ", output)
        
        if binary_output is not None:
            binary_mask = np.where(self.masks !=0, 1, 0).astype(np.uint8)
            tifffile.imwrite(binary_output, binary_mask)
            print("Saved binary mask to: ", binary_output)

    
    def plot_qc(self, qc_output):
        if qc_output is not None:
            fig = plt.figure(figsize=(12,5))
            plot.show_segmentation(fig, self.dapi_channel, self.masks, self.flows[0])
            plt.tight_layout()
            plt.savefig(qc_output)
            plt.close()
            print("Saved QC png to: ", qc_output)

        
        

def parse_args():
    parser = argparse.ArgumentParser(description='Find Nuclei using cellpose on DAPI channel',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-i', '--input', 
                       required=True,
                       help='input image filepath')
    
    parser.add_argument('-o', '--output',
                       required=True, 
                       help='output mask filepath')
    
    parser.add_argument('-b', '--binary_output',
                       required=False,
                       default=None,
                       help='output binary mask filepath if requested')
    
    parser.add_argument('-q', '--qc_output',
                       required=False,
                       default=None,
                       help='QC Output PNG')
    
    parser.add_argument('-n', '--dapi_index',
                       required=False,
                       default=0,
                       type=int,
                       help='Index of DAPI in image (zero-based indexing)')
    
    parser.add_argument('-c', '--checkpoint',
                       required=False,
                       default="/group/kidn3/Emma/MODEL_CHECKPOINTS/sam2.1_hiera_large.pt",
                       help='Model checkpoint')
    
    parser.add_argument('-m', '--model_config',
                       required=False,
                       default="configs/sam2.1/sam2.1_hiera_l.yaml",
                       help='Model configuration yaml')
    
    parser.add_argument('-d', '--device',
                       required=False,
                       default="cuda",
                       choices=["cuda", "cpu"],
                       help='Device')
    
    parser.add_argument('-p', '--diameter',
                       required=False,
                       default=30,
                       type=int,
                       help='Diameter for nuclei segmentation')
    
    return parser.parse_args()


def main():
    args = parse_args()
    segmentor = NucleiSegmentor(args.input, dapi_index=args.dapi_index)
    segmentor.find_nuclei_cellpose(args.diameter)
    segmentor.save_outputs(args.output, args.binary_output)
    segmentor.plot_qc(args.qc_output)

if __name__ == "__main__":
    main()

