from time import time

import matplotlib.pyplot as plt
import scipy.optimize as spo

from pyrt.data.data_trots import *
from pyrt.tools import *
from pyrt.optimization.model import *

import gurobipy as grb







class vmat_mip(model_base):
    def __init__(self, input_dict,modality='imrt', run_title='default'):
        super(self.__class__,self).__init__(input_dict,modality=modality)
        self.run_title = run_title
        self.model_params = input_dict['model_params'].copy()

        self.m = grb.Model()

        self.build_model()



    def build_model(self):
        # generate optimization metadata you need, can break into functions like you did before (building vars, dose vars, aper vars, etc)
        # you have self.data as the data object


        pass



    def build_dose_constraints(self):
        pass


    def optimize(self): # probably need some inputs/default params like in IMRT one


        # write optimization code here (probably self.m.optimize(), then some data extraction)


        # save apertures (or the indices of apertures) of solution

        self.dose_dict[self.run_title] = [np.copy(self.current_dose_per_structure[s]) for s in range(len(self.data.structures))]
        self.obj_dict[self.run_title] = self.current_obj


    def calc_dose_from_variables(self):
        # given apertures, calculate dose
        #  will be useful in long run

        pass




    def calc_obj_grad(self, beamlet_intensities):

        pass




