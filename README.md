# hawk_to_pypulseq
A conversion script for RTHAWK pulse sequences (apd + spv files) into the pypulseq format for data sharing.

Tested and built on python 3.8.10.

# Instructions
## Create a virtual environment
1. `python3 -m venv venv`
2. `source venv/bin/activate`

## Install required packages using pip
`pip install -r requirements.txt`

## Run code:
`python3 run_hawk_to_pulseq_convertor.py`

## Dynamic Imaging Reconstruction
Submodule for a MATLAB ISMRMRD reconstruction using FIRE is linked to this repository for it to be a main resource for ISMRM abstract 2402: Open-source dynamic MRI workflow for reproducible research
