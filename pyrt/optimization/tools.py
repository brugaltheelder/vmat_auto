
import numpy as np
from pyrt.data.tools import control_point
from pyrt.data.data_trots import patient_data

def save_weights_from_input_dict(model, target_weights_label='target_weights', oar_weights_label='oar_weights'):
    # weights
    model.weights_per_structure = [np.ones(s.num_vox) for s in model.data.structures]
    for s_ind in range(len(model.data.structures)):
        if model.data.structures[s_ind].is_target:
            if model.data.structures[s_ind].name in model.model_params[target_weights_label].keys():
                model.weights_per_structure[s_ind] *= model.model_params[target_weights_label][
                    model.data.structures[s_ind].name]
            else:
                model.weights_per_structure[s_ind] *= model.model_params[target_weights_label]['default']
        else:
            if model.data.structures[s_ind].name in model.model_params[oar_weights_label].keys():
                model.weights_per_structure[s_ind] *= model.model_params[oar_weights_label][
                    model.data.structures[s_ind].name]
            else:
                model.weights_per_structure[s_ind] *= model.model_params[oar_weights_label]['default']


class aperture(object):
    def __init__(self, data, CP, aper_left_pos = [], aper_right_pos = [], aper_intensity = 1. , set_open_aper=False):
        ### aper shape details
        assert(isinstance(CP,control_point))
        assert (isinstance(data, patient_data))
        if set_open_aper:
            self.left_leaves = [CP.left_leaf_position[r] for r in range(CP.num_rows)]
            self.right_leaves = [CP.left_leaf_position[r] + CP.width_per_row[r] for r in range(CP.num_rows)]
            self.intensity = 1.
        else:
            self.left_leaves = aper_left_pos[:]
            self.right_leaves = aper_right_pos[:]
            self.intensity = aper_intensity

        self.build_Dkj(CP,data)



    def build_Dkj(self,CP,data):
        assert (isinstance(CP, control_point))
        assert (isinstance(data, patient_data))
        for s in range(len(data.structures)):
            dense_dkj = np.zeros(data.structures[s].num_vox)
            for r in range(CP.num_rows):
                for i in range(self.left_leaves[r], self.right_leaves[r]):
                    beamlet_index = CP.left_leaf_index[i]
            self.Dkj_per_structure = None


    def calc_dose(self,x=None):
        pass