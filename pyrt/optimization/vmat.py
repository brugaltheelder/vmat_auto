from time import time

import matplotlib.pyplot as plt
import scipy.optimize as spo

from pyrt.data.data_trots import *
from pyrt.tools import *
from pyrt.optimization.model import *
from pyrt.optimization.tools import *

import gurobipy as grb






class conformal_arc(model_base):

    def __init__(self, input_dict,modality='conf_arc', run_title='default'):
        super(self.__class__,self).__init__(input_dict,modality=modality)
        self.run_title = run_title
        self.model_params = input_dict['model_params'].copy()

        self.num_apers = self.data.num_control_points

        self.apertures = []
        self.build_model()
        self.generate_apertures()




    def build_model(self, target_weights_label='target_weights', oar_weights_label='oar_weights'):
        # variables
        self.aper_intensities = np.zeros(self.num_apers)
        self.grad = np.zeros(self.num_apers)
        self.zhat = [np.zeros(s.num_vox) for s in self.data.structures]

        # parameters
        # thresholds
        self.thresholds_per_structure = [np.ones(s.num_vox)*s.rx for s in self.data.structures]
        # weights
        self.weights_per_structure = [np.ones(s.num_vox) for s in self.data.structures]
        save_weights_from_input_dict(self)



    def generate_apertures(self):
        for a in range(self.num_apers):
            self.apertures.append(aperture(self.data, self.data.control_points[a],set_open_aper=True))



    def optimize(self,startingX=None, UB=None, display=False, ftol=1e-5, gtol=1e-5):
        if startingX is not None:
            x0 = np.copy(startingX)
        else:
            x0 = np.zeros(self.num_apers)

        print 'Solver Starting...'
        start = time()
        res = spo.minimize(self.calc_obj_grad, x0=x0, method='L-BFGS-B', jac=True, bounds=[(0, UB) for i in
                        xrange(self.num_apers)], options={'ftol': ftol, 'gtol': gtol, 'disp': display})

        print '{} model solved in {} seconds'.format(self.modality,time() - start)

        self.aper_intensities, self.current_obj =res['x'], res['fun']
        for a in range(len(self.apertures)):
            self.apertures[a].intensity = self.aper_intensities[a]
        self.calc_dose_from_variables(self.aper_intensities)


        self.dose_dict[self.run_title] = [np.copy(self.current_dose_per_structure[s]) for s in range(len(self.data.structures))]
        self.obj_dict[self.run_title] = self.current_obj


    def calc_dose_from_variables(self,x=None):
        self.current_dose_per_structure = [np.zeros(self.data.structures[s].num_vox) for s in range(len(self.data.structures))]

        if x is not None:
            for a in range(len(self.apertures)):
                for s in range(len(self.data.structures)):
                    self.current_dose_per_structure[s] += self.apertures[a].Dkj_per_structure[s] * x[a]
        else:
            for a in range(len(self.apertures)):
                for s in range(len(self.data.structures)):
                    self.current_dose_per_structure[s] += self.apertures[a].Dkj_per_structure[s] * self.apertures[s].intensity


    def get_obj_update_zhat(self):
        obj = 0.
        for s in range(len(self.data.structures)):
            dose_deviation = np.array(self.current_dose_per_structure[s]-self.thresholds_per_structure[s])
            obj += float(self.weights_per_structure[s].dot(dose_deviation**2))

            self.zhat[s] = np.multiply(self.weights_per_structure[s], dose_deviation)
        return obj


    def calc_obj_grad(self, aper_intensities):
        self.calc_dose_from_variables(aper_intensities)

        self.grad.fill(0.)
        for s in range(len(self.data.structures)):
            self.zhat[s].fill(0.)

        obj = 0.0
        obj += self.get_obj_update_zhat()
        for a in range(len(self.apertures)):
            for s in range(len(self.data.structures)):
                self.grad[a] += self.apertures[a].Dkj_per_structure[s].dot(self.zhat[s])
                # self.grad += self.data.structures[s].Dij.dot(self.zhat[s])

        return obj, self.grad




class vmat_mip(model_base):
    def __init__(self, input_dict,modality='imrt', run_title='default',build_model = True):
        super(self.__class__,self).__init__(input_dict,modality=modality)
        self.run_title = run_title
        self.model_params = input_dict['model_params'].copy()
        self.weights_per_structure = None

        self.m = grb.Model()

        save_weights_from_input_dict(self)

        self.apertures_per_cp = [[] for cp in range(self.data.num_control_points)]

        self.build_apertures()


        if build_model:
            self.build_model()

    def build_apertures(self, aper_types_list=['conf']):

        if 'conf' in aper_types_list:
            # execute conf aper building
            for cp in range(self.data.num_control_points):

                self.apertures_per_cp[cp].append(aperture(self.data, self.data.control_points[cp], set_open_aper=True))


    def build_model(self):
        # generate optimization metadata you need, can break into functions like you did before (building vars, dose vars, aper vars, etc)
        # you have self.data as the data object
        print_in_box('Building Gurobi Model')
        print 'Initializing dose and aperture variables'
        self.generate_dose_variables()
        self.generate_aper_variables()

        self.generate_thresholds()
        print 'Initializing objective variables'
        self.generate_objective_variables()
        self.generate_aperture_per_cp()

        # print 'Building aperture network'
        # self.generate_network()

        print 'Creating constraints'
        self.generate_bilinear_constraints()
        self.build_dose_constraints()
        self.generate_objective_constraints()
        self.generate_objective()

        print_in_box('Gurobi Model Constructed')


    def generate_dose_variables(self):
        self.dose_var = [self.m.addVars(self.data.structures[s].num_vox, name='z_{}_'.format(self.data.structures[s].name.replace(' ','_'))) for s in range(len(self.data.structures))]
        self.m.update()

    def generate_aper_variables(self):
        self.aper_intensity_var = [[self.m.addVar(lb=0,name='y_{}_{}'.format(cp,a)) for a in range(len(self.apertures_per_cp[cp]))] for cp in range(self.data.num_control_points)]
        self.aper_binary_var = [[self.m.addVar(lb=0, vtype=grb.GRB.BINARY,name='x_{}_{}'.format(cp,a)) for a in range(len(self.apertures_per_cp[cp]))] for cp in range(self.data.num_control_points)]
        self.m.update()

    def generate_thresholds(self):
        self.thresholds = [[self.data.structures[s].rx for v in range(self.data.structures[s].num_vox)] for s in range(len(self.data.structures))]
        self.m.update()

    def generate_objective_variables(self):
        self.obj_var = [self.m.addVars(self.data.structures[s].num_vox, name='h_{}_'.format(self.data.structures[s].name.replace(' ','_'))) for s in range(len(self.data.structures))]
        self.m.update()

    def build_dose_constraints(self):
        for s in range(len(self.data.structures)):
            print 'Building dose constraint for {}'.format(self.data.structures[s].name)
            for v in range(self.data.structures[s].num_vox):
                lin1 = grb.LinExpr()
                for cp in range(self.data.num_control_points):
                    for a in range(len(self.apertures_per_cp[cp])):
                        self.apertures_per_cp[cp][a].Dkj_per_structure[s][v]
                    lin1 += self.apertures_per_cp[cp][a].Dkj_per_structure[s][v] * self.aper_intensity_var[cp][a]
                self.m.addConstr(self.dose_var[s][v], grb.GRB.EQUAL, lin1, name='Dose_Constraint_{}'.format(self.data.structures[s].name.replace(' ','_')))

        self.m.update()

    # Select only one aperture per beam
    def generate_aperture_per_cp(self):
        for cp in range(self.data.num_control_points):
            self.m.addConstr(sum(self.aper_binary_var[cp][a] for a in range(len(self.apertures_per_cp[cp]))), grb.GRB.EQUAL, self.model_params['aper_limit'],name='One_Aperture_Per_Beam_{}'.format(cp))
        self.m.update()

    def generate_bilinear_constraints(self):
        for cp in range(self.data.num_control_points):
            for a in range(len(self.apertures_per_cp[cp])):
                self.m.addConstr(self.aper_intensity_var[cp][a], grb.GRB.LESS_EQUAL,self.aper_binary_var[cp][a] * self.model_params['max_intensity'], name='Bilinear_Upper_Bound_[{}][{}]'.format(cp, a))
                self.m.addConstr(self.aper_intensity_var[cp][a], grb.GRB.GREATER_EQUAL, self.aper_binary_var[cp][a] * self.model_params['min_intensity'], name='Bilinear_Lower_Bound_[{}][{}]'.format(cp, a))
        self.m.update()

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


        # write optimization code here, then some data extraction)
        self.m.optimize()

        # save apertures (or the indices of apertures) of solution

        for s in range(len(self.data.structures)):
            for v in range(self.data.structures[s].num_vox):
                self.current_dose_per_structure[s][v] = self.dose_var[s][v].x

        # todo extract other variable info (aperture intensity)

        self.dose_dict[self.run_title] = [np.copy(self.current_dose_per_structure[s]) for s in range(len(self.data.structures))]
        self.obj_dict[self.run_title] = self.current_obj


    def calc_dose_from_variables(self):
        # given apertures, calculate dose
        #  will be useful in long run

        pass




    def calc_obj_grad(self, beamlet_intensities):

        pass




