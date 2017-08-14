from time import time

import matplotlib.pyplot as plt
import scipy.optimize as spo

from pyrt.data.data_trots import *
from pyrt.tools import *
from pyrt.optimization.model import *
from pyrt.optimization.tools import *

import gurobipy as grb







class vmat_mip(model_base):
    def __init__(self, input_dict,modality='imrt', run_title='default'):
        super(self.__class__,self).__init__(input_dict,modality=modality)
        self.run_title = run_title
        self.model_params = input_dict['model_params'].copy()
        self.weights_per_structure = None

        self.m = grb.Model()

        save_weights_from_input_dict(self)


        self.build_model()



    def build_model(self):
        # generate optimization metadata you need, can break into functions like you did before (building vars, dose vars, aper vars, etc)
        # you have self.data as the data object

        self.generate_dose_variables()
        # self.generate_aper_variables()
        self.generate_thresholds()
        self.generate_objective_variables()
        # self.generate_one_aperture()
        # self.generate_network()
        # self.generate_bilinear_constraints()
        # self.generate_dose_constraints()
        self.generate_objective_constraints()
        self.generate_objective()

        pass

    def generate_dose_variables(self):
        self.dose_var = [self.m.addVars(self.data.structures[s].num_vox, name='z_{}_'.format(self.data.structures[s].name.replace(' ','_'))) for s in range(len(self.data.structures))]
        self.m.update()

    def generate_thresholds(self):
        self.thresholds = [[self.data.structures[s].rx for v in range(self.data.structures[s].num_vox)] for s in range(len(self.data.structures))]
        self.m.update()

    def generate_objective_variables(self):
        self.obj_var = [self.m.addVars(self.data.structures[s].num_vox, name='h_{}_'.format(self.data.structures[s].name.replace(' ','_'))) for s in range(len(self.data.structures))]
        self.m.update()

    def build_dose_constraints(self):
        # for s in range(len(self.data.structures)):
        #     for v in range(self.data.structures[s].num_vox):
        #         lin1 = grb.LinExpr()
        #         for cp in range(model.data.num_control_points):
        #         for aper_index in range(self.data.apertures_per_beam_dict[beam_index].num_apers):
        #             lin1 += self.data.apertures_per_beam_dict[beam_index].Dbaj[aper_index, voxel] * self.data.apertures_per_beam_dict[beam_index].intensity_variables[aper_index]
        #     self.m.addConstr(self.dose_var[voxel], grb.GRB.EQUAL, lin1, name='Dose_Constraint_{}'.format(voxel))
        #
        # self.m.update()
        pass

    # h >= z - Rx, h >=0, and h >= Rx -z if tumor
    def generate_objective_constraints(self):
        for s in range(len(self.data.structures)):
            for v in range(self.data.structures[s].num_vox):
                self.m.addConstr(self.obj_var[s][v], grb.GRB.GREATER_EQUAL, 0.)
                self.m.addConstr(self.obj_var[s][v], grb.GRB.GREATER_EQUAL, self.dose_var[s][v] - self.thresholds[s][v])
                if self.data.structures[s].is_target:
                    self.m.addConstr(self.obj_var[s][v], grb.GRB.GREATER_EQUAL, self.thresholds[s][v] - self.dose_var[s][v])

        self.m.update()

    def generate_objective(self):
        lin2 = grb.LinExpr()
        for s in range(len(self.data.structures)):
            for v in range(self.data.structures[s].num_vox):
                lin2 += self.obj_var[s][v] * self.weights_per_structure[s][v]

        self.m.setObjective(lin2, grb.GRB.MINIMIZE)
        self.m.update()

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




