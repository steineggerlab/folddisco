# Read pickle tensor file and convert to safetensor file
# Usage: python3 convert_pt_to_safetensor.py <pickle_file> <safetensor_file>
#.       python3 convert_pt_to_safetensor.py <pickle_file>

import sys
import torch
from safetensors.torch import save_file

def main():
    # Get input and output file names
    input_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    else:
        output_file = input_file + ".safetensors"

    # Read pickle file
    with open(input_file, 'rb') as f:
        tensor = torch.load(f, map_location=torch.device('cpu'))
        # Save tensor to safetensor file
        save_file(tensor.state_dict(), output_file)

if __name__ == "__main__":
    main()