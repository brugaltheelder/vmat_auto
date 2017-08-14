from pyrt.data.data_trots import *


# Find work directory
cwd = '/Users/jtmargo/Desktop/Josh/Clemson University/Research/Treatment Plan/TROTS Data/Prostate_VMAT/'

# Load patient information
filename = 'Prostate_VMAT_101.mat'
# f = h5py.File(cwd + filename, 'r')

input_dict = {
    'cwd': cwd,
    'filename': filename,
    'Rx': {'PTV': 79.56, 'PTV Vesicles': 72.2},
}

model = patient_data(input_dict)
# print model.structure_names