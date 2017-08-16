import numpy as np
from pyrt.data.tools import control_point_vmat
from pyrt.data.data_trots import patient_data
from pyrt.tools import *


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


def is_connected(aperture1, aperture2, max_distance_per_cp):
    if not aperture1.num_rows == aperture2.num_rows:
        print_in_box('NUMBER OF ROWS BETWEEN {} and {} DOES NOT MATCH'.format(aperture1.cp_number, aperture2.cp_number))
    else:
        total_rows = aperture1.num_rows
    var = abs(aperture1.cp_number - aperture2.cp_number) * max_distance_per_cp
    if all([abs(aperture1.left_leaf_position[r] - aperture2.left_leaf_position[r]) <= var for r in
            range(total_rows)]) and all(
            [abs(aperture1.right_leaf_position[r] - aperture1.right_leaf_position[r]) <= var for r in
             range(total_rows)]):
        return True
    else:
        return False


def aper_gen_given_dose(dose_per_struct, data, CP):
    # todo troy write this function

    return None


class aperture(object):
    def __init__(self, data, CP, aper_left_pos=[], aper_right_pos=[], aper_intensity=1., set_open_aper=False,
                 beamlet_override=None):
        ### aper shape details
        assert (isinstance(CP, control_point_vmat))
        assert (isinstance(data, patient_data))
        self.cp_number = CP.cp_number
        self.num_rows = CP.num_rows
        if set_open_aper:

            self.left_leaf_position = [CP.left_leaf_position[r] for r in range(CP.num_rows)]
            self.right_leaf_position = [CP.left_leaf_position[r] + CP.width_per_row[r] for r in range(CP.num_rows)]
            self.intensity = 1.
        elif beamlet_override is not None:
            # todo troy beamlet override
            pass


        else:
            self.left_leaf_position = aper_left_pos[:]
            self.right_leaf_position = aper_right_pos[:]
            self.intensity = aper_intensity

        self.build_Dkj(CP, data)

    def build_Dkj(self, CP, data):
        assert (isinstance(CP, control_point_vmat))
        assert (isinstance(data, patient_data))

        self.Dkj_per_structure = [np.zeros(s.num_vox) for s in data.structures]
        self.beamlet_members = []

        for r in range(CP.num_rows):
            for i in range(self.left_leaf_position[r], self.right_leaf_position[r]):
                self.beamlet_members.append(
                    CP.initial_beamlet_index + CP.left_leaf_index[r] + (i - CP.left_leaf_position[r]))

        for s in range(len(data.structures)):
            self.Dkj_per_structure[s] += np.asarray(
                data.structures[s].Dij[np.array(self.beamlet_members)].sum(axis=0)).flatten()

    def get_fluence_map(self):

        fluence = np.zeros(len(self.left_leaf_position), max(self.right_leaf_position))
        for r in range(len(self.left_leaf_position)):
            fluence[self.left_leaf_position[r]:self.right_leaf_position[r]] = self.intensity
        return fluence

    def calc_dose(self, x=None):
        pass
