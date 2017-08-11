# this file will have "model" structured files

import numpy as np
import scipy.optimize as spo
import scipy.sparse as sps

from data import *




class model(object):
    def __init__(self, input_dict, modality='undefined'):
        self.read_in_data(input_dict)
        self.structure_dose = [np.zeros(s.num_vox) for s in self.data.structures]


    def read_in_data(self,input_dict):
        self.data = patient_data(input_dict)

    def generate_optimization_variables(self):
        print 'Implement in child classes'
        pass

    def optimize(self):
        print 'Implement in child classes'
        pass

    def calc_dose_from_variables(self):
        print 'Implement in child classes'
        pass

    def plot_dvh(self):
        print 'Implement in child classes'
        pass




class imrt_fmo(model):
    def __init__(self, input_dict,modality='imrt'):
        super(self.__class__,self).__init__(input_dict,modality=modality)



    def generate_optimization_variables(self):
        # variables
        self.beamlet_intensities = np.zeros(self.data.num_beamlets)

        # parameters
        # thresholds

        # weights

    def optimize(self):
        pass

    def calc_obj_grad(self):
        pass






