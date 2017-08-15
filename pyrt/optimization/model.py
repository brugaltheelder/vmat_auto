# this file will have "model" structured files

from time import time

import matplotlib.pyplot as plt
import scipy.optimize as spo

from pyrt.data.data_trots import *
from pyrt.tools import plot_DVH, print_in_box


class model_base(object):
    def __init__(self, input_dict, modality='undefined'):
        self.modality = modality
        print_in_box('Reading in data')
        self.read_in_data(input_dict, modality)
        print_in_box('Data reading completed')
        self.current_dose_per_structure = [np.zeros(s.num_vox) for s in self.data.structures]
        self.current_obj = 0.
        self.dose_dict = {}
        self.obj_dict = {}

    def read_in_data(self,input_dict, modality):
        self.data = patient_data(input_dict, modality)

    def build_model(self):
        print 'Implement in child classes'
        pass

    def optimize(self):
        print 'Implement in child classes'
        pass

    def calc_dose_from_variables(self,x=None):
        print 'Implement in child classes'
        pass

    def save_current_dose(self,run_tag):
        self.dose_dict[run_tag] = [np.copy(self.current_dose_per_structure[s]) for s in range(len(self.data.structures))]


    def plot_DVH(self, saveName='', showPlots=False, saveDVH=False, run_tag=None,num_bins = 100):
        plot_DVH(self, saveName=saveName,showPlots=showPlots,saveDVH=saveDVH, run_tag=run_tag,num_bins=500)


    def serialize_dose(self):
        pass















