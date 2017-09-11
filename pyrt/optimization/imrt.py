
from time import time

import matplotlib.pyplot as plt
import scipy.optimize as spo

from pyrt.data.data_trots import *
from pyrt.tools import *
from pyrt.optimization.model import *
from pyrt.optimization.tools import *


class imrt_fmo(model_base):
    def __init__(self, input_dict,modality='imrt', run_title='default', tora = False, tora_dose=None):
        super(self.__class__,self).__init__(input_dict,modality=modality)
        self.run_title = run_title
        self.model_params = input_dict['model_params'].copy()
        self.tora = tora
        if self.tora:
            if tora_dose is None:
                print 'SET TORA DOSE!'
                exit()
            else:
                #todo write dose sorting function
                self.tora_dose = [d.copy() for d in tora_dose]
                self.sorted_target_dose_dict =  [np.ones(s.num_vox) for s in self.data.structures]
        self.build_model()



    def build_model(self, target_weights_label='target_weights', oar_weights_label='oar_weights'):
        # variables
        self.beamlet_intensities = np.zeros(self.data.num_beamlets)
        self.grad = np.zeros(self.data.num_beamlets)
        self.zhat = [np.zeros(s.num_vox) for s in self.data.structures]

        # parameters
        # thresholds
        self.thresholds_per_structure = [np.ones(s.num_vox)*s.rx for s in self.data.structures]
        if self.tora:
            pass


        # weights
        self.weights_per_structure = None
        if self.tora:
            self.gen_weights_tora()
        else:
            save_weights_from_input_dict(self, target_weights_label=target_weights_label,oar_weights_label=oar_weights_label)


    def gen_tora_threshold(self, scalar = 0., sorted_index_dict = None):
        self.threshold_primary_over = [np.ones(s.num_vox) for s in self.data.structures]
        self.threshold_primary_under = [np.ones(s.num_vox) for s in self.data.structures]

        # resort target_dose_dist
        if sorted_index_dict is not None:
            for s in range(self.data.structures):
                if self.data.structures[s].is_target:
                    self.threshold_primary_over[s][sorted_index_dict[s]] = np.maximum((1 - scalar) * self.sorted_target_dose_dict[
                        s] + scalar * self.thresholds_per_structure[s][sorted_index_dict[s]],self.thresholds_per_structure[s][sorted_index_dict[s]])
                    self.threshold_primary_under[s][sorted_index_dict[s]] = np.minimum((1 - scalar) * self.sorted_target_dose_dict[
                        s] + scalar * self.thresholds_per_structure[s][sorted_index_dict[s]],self.thresholds_per_structure[s][sorted_index_dict[s]])
                else:
                    self.threshold_primary_over[s][sorted_index_dict[s]] = self.sorted_target_dose_dict[s]
                    self.threshold_primary_under[s][sorted_index_dict[s]] = self.sorted_target_dose_dict[s]

        else:
            for s in range(self.data.structures):
                self.threshold_primary_over[s] = (1-scalar) * self.tora_dose[s] + scalar * self.thresholds_per_structure[s]
                self.threshold_primary_under[s] = (1-scalar) * self.tora_dose[s] + scalar * self.thresholds_per_structure[s]






    def gen_weights_tora(self, epsilon=0.001):
        self.weights_per_structure = [np.ones(s.num_vox) for s in self.data.structures]
        self.primary_weights_over  = [np.ones(s.num_vox) for s in self.data.structures]
        self.primary_weights_under = [np.ones(s.num_vox) for s in self.data.structures]
        for s_ind in range(len(self.data.structures)):
            if self.data.structures[s_ind].is_target:
                self.primary_weights_under[s_ind] *= 1./self.data.structures[s_ind].num_vox
            else:
                self.primary_weights_under[s_ind] *= 0.
            self.primary_weights_over *= 1./self.data.structures[s_ind].num_vox

            # meanStructureDoseDeviation = max(np.average((self.ideal_dose_dist[self.structure_index_dict[s]] -
            #                                              self.target_dose_dist[self.structure_index_dict[s]]) ** 2), 1.)
            #




    def optimize(self,startingX=None, UB=None, display=False, ftol=1e-5, gtol=1e-5, tora=False):
        if startingX is not None:
            x0 = np.copy(startingX)
        else:
            x0 = np.zeros(self.data.num_beamlets)

        if tora==False:
            start = time()
            res = spo.minimize(self.calc_obj_grad, x0=x0, method='L-BFGS-B', jac=True, bounds=[(0, UB) for i in
                            xrange(self.data.num_beamlets)], options={'ftol': ftol, 'gtol': gtol, 'disp': display})

            print '{} model solved in {} seconds'.format(self.modality,time() - start)

            self.beamlet_intensities, self.current_obj =res['x'], res['fun']
            self.calc_dose_from_variables()

            self.dose_dict[self.run_title] = [np.copy(self.current_dose_per_structure[s]) for s in range(len(self.data.structures))]
            self.obj_dict[self.run_title] = self.current_obj
        else:
            # tora stuff
            pass


    def calc_dose_from_variables(self,x=None):
        if x is not None:
            for s in range(len(self.data.structures)):
                self.current_dose_per_structure[s] = self.data.structures[s].Dij.transpose().dot(x)

        else:
            for s in range(len(self.data.structures)):
                self.current_dose_per_structure[s] = self.data.structures[s].Dij.transpose().dot(self.beamlet_intensities)


    def get_obj_update_zhat(self):
        obj = 0.
        for s in range(len(self.data.structures)):
            dose_deviation = np.array(self.current_dose_per_structure[s]-self.thresholds_per_structure[s])
            obj += float(self.weights_per_structure[s].dot(dose_deviation**2))
            self.zhat[s] = np.multiply(self.weights_per_structure[s], dose_deviation)
        return obj


    def calc_obj_grad(self, beamlet_intensities):
        self.calc_dose_from_variables(beamlet_intensities)

        self.grad.fill(0.)
        for s in range(len(self.data.structures)):
            self.zhat[s].fill(0.)

        obj = 0.0
        obj += self.get_obj_update_zhat()
        for s in range(len(self.data.structures)):
            self.grad += self.data.structures[s].Dij.dot(self.zhat[s])

        return obj, self.grad

