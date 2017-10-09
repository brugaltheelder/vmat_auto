from pyrt.optimization.vmat import *
from pyrt.research.run_functions import *

# Josh Inputs - Prostate
# Find work directory
# cwd = '/Users/troy/Dropbox/CAP Group/TROTS/Prostate_VMAT/'
working_directory = 'outputs_josh/'
path1 = '/Users/jtmargo/Desktop/Josh/Clemson University/Research/Treatment Plan/TROTS Data/Prostate_VMAT/'
path2 = '/Users/jtmargo/Desktop/Josh/Clemson University/Research/Treatment Plan/TROTS Data/Head-and-Neck/'


cwd = path1

if cwd == path1:
    back_proj_list_of_dicts = [
       {'PTV':30.,'Bladder':15.,'default':0., 'threshold':0.},
        {'PTV':30.,'Rectum':15.,'default':0., 'threshold':0.},
    #     {'PTV':5,'default':1, 'threshold':0.},
        {'PTV':100.,'default':1., 'threshold':0.},
        {'PTV':1.,'default':0., 'threshold':0.}
    ]
    vmat_model_params = {
        'target_weights':{'PTV':1000., 'default':500.},
        'oar_weights':{'Rectum':10., 'Bladder':10., 'default':1.},
        'max_intensity':1000.,
        'min_intensity': 1.,
        'aper_limit': 1.,
        'max_distance_per_cp': 100,
        'cp_redundancy': 1,
        'back_projection_dicts':back_proj_list_of_dicts
    }

    input_dict = {
        'cwd': cwd,
        'figure_directory':working_directory,
        'aper_types_list': ['back_proj'],
        'filename': None,
        'Rx': {'PTV': 79.56, 'PTV Vesicles': 72.2},
        'model_params':vmat_model_params
    }

if cwd == path2:
    back_proj_list_of_dicts = [
       {'PTV':30.,'PTV 0-46 Gy':30.,'Brainstem':15.,'default':0., 'threshold':0.},
        {'PTV':30.,'PTV 0-46 Gy':30.,'Spinal Cord':15.,'default':0., 'threshold':0.},
    #     {'PTV':5,'default':1, 'threshold':0.},
        {'PTV':100.,'PTV 0-46 Gy':100.,'default':1., 'threshold':0.},
        {'PTV':1.,'PTV 0-46 Gy':1.,'default':0., 'threshold':0.}
    ]

    vmat_model_params = {
        'target_weights':{'PTV':1000., 'PTV 0-46 Gy':1000., 'default':500.},
        'oar_weights':{'Brainstem':10., 'Spinal Cord':10., 'default':1.},
        'max_intensity':1000.,
        'min_intensity': 1.,
        'aper_limit': 1.,
        'max_distance_per_cp': 100,
        'cp_redundancy': 1,
        'back_projection_dicts':back_proj_list_of_dicts
    }

    input_dict = {
        'cwd': cwd,
        'figure_directory':working_directory,
        'aper_types_list': ['back_proj'],
        'filename': None,
        'Rx': {'PTV':47.15,'PTV 0-46 Gy': 47.15},
        'model_params':vmat_model_params
    }


#Run individual case

# Load patient information
input_dict['filename'] = 'Prostate_VMAT_210.mat'
# input_dict['filename'] = 'Head-and-Neck_05.mat'
# input_dict['filename'] = 'HN01_FullVMAT.mat'
input_dict['case_directory'] = input_dict['filename'][0:-4]+'/'


run_case(input_dict)