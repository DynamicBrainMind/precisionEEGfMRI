import pickle
from scipy.io import savemat
import glob
import os

# Folder containing the .pkl files
folder = '/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/formatted/full'

# Find all .pkl files
files = glob.glob(os.path.join(folder, '*.pkl'))

for pkl_file in files:

    # Load pickle file
    with open(pkl_file, 'rb') as f:
        data = pickle.load(f)

    # Create output filename
    base = os.path.splitext(pkl_file)[0]
    mat_file = base + '.mat'

    # Save as MATLAB file
    savemat(mat_file, {'data': data})

    print(f'Saved: {mat_file}')