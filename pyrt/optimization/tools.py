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


def gen_weighted_mask(dose_per_struct, data):
    target_mask = [np.ones(data.structures[s].num_vox) * dose_per_struct['default'] for s in
                   range(len(data.structures))]
    for s in range(len(data.structures)):
        if data.structures[s].name in dose_per_struct.keys():
            factor = 1.
            if data.structures[s].is_target:
                factor = -1.

            target_mask[s] = factor * dose_per_struct[data.structures[s].name] * np.ones(data.structures[s].num_vox) / \
                             data.structures[s].num_vox
        elif data.structures[s].is_target:
            target_mask[s] = -1. * target_mask[s]* np.ones(data.structures[s].num_vox) / \
                             data.structures[s].num_vox
        else:
            target_mask[s] = 1. * target_mask[s] * np.ones(data.structures[s].num_vox) / \
                             data.structures[s].num_vox
    return target_mask


def aper_gen_given_dose(dose_per_struct, data, CP, mask=None):
    if mask is not None:
        target_mask = mask
    else:
        target_mask = gen_weighted_mask(dose_per_struct, data)

    beamlet_usefulness = np.zeros(CP.number_beamlets)

    for s in range(len(data.structures)):
        beamlet_usefulness += data.structures[s].Dij[np.arange(CP.initial_beamlet_index, CP.final_beamlet_index)].dot(
            target_mask[s])

    value, beamlets = pricing_problem_beam_aper(CP, beamlet_usefulness)

    return aperture(data, CP, beamlet_override=beamlets),beamlet_usefulness




def pricing_problem_beam_aper(CP, beamlet_gradient, threshold=0.):
    beamlets = []

    # todo implement thresholding?

    worth = 0
    assert (isinstance(CP, control_point_vmat))
    # for each row
    for r in range(CP.num_rows):
        maxSoFar, maxEndingHere, lE, rE, = 0, 0, 0, 0
        for i in range(CP.width_per_row[r]):
            maxEndingHere += beamlet_gradient[i + CP.left_leaf_index[r]]
            if maxEndingHere > 0:
                maxEndingHere, lE, rE = 0, i + 1, i + 1
            if maxSoFar > maxEndingHere:
                maxSoFar, rE = maxEndingHere, i + 1

        for i in range(lE, rE):
            beamlets.append(i + CP.left_leaf_index[r])

        worth += maxSoFar
    return worth, beamlets


class aperture(object):
    def __init__(self, data, CP, aper_left_pos=[], aper_right_pos=[], aper_intensity=1., set_open_aper=False,
                 beamlet_override=None):
        ### aper shape details
        assert (isinstance(CP, control_point_vmat))

        assert (isinstance(data, patient_data))
        self.cp_number = CP.cp_number
        self.org_cp_number = CP.org_cp_number
        self.num_rows = CP.num_rows

        if set_open_aper:
            self.left_leaf_position = [CP.left_leaf_position[r] for r in range(CP.num_rows)]
            self.right_leaf_position = [CP.left_leaf_position[r] + CP.width_per_row[r] for r in range(CP.num_rows)]
            self.intensity = 1.
        elif beamlet_override is not None:
            dummy_field = np.zeros_like(CP.field)
            beamlet_counter = 1
            for b in beamlet_override:
                dummy_field[CP.field_position[int(b)][0], CP.field_position[int(b)][1]] = beamlet_counter
                beamlet_counter += 1
            self.build_leaf_metadata(dummy_field, CP)
            self.intensity = aper_intensity



        # elif beamlet_override is not None and False:
        #
        #
        #     row_counter = 0
        #     self.left_leaf_position = [int((CP.left_leaf_position[r]+CP.width_per_row[r])/2 ) for r in range(CP.num_rows)]
        #     self.right_leaf_position = [int((CP.left_leaf_position[r]+CP.width_per_row[r])/2 ) for r in range(CP.num_rows)]
        #     self.intensity = aper_intensity
        #     current_row = int(CP.field_position[beamlet_override[0]][0])
        #     self.left_leaf_position[0] = int(CP.field_position[beamlet_override[0]][1])
        #     print current_row
        #     print self.left_leaf_position
        #     print self.right_leaf_position
        #     total_beamlets = 0
        #     break_bool = False
        #     for b in range(1,len(beamlet_override)):
        #
        #         x,y = tuple(CP.field_position[beamlet_override[b]])
        #
        #         if int(x)>current_row:
        #             self.right_leaf_position[row_counter] = int(CP.field_position[beamlet_override[b - 1]][1]) + 1
        #             row_counter+=1
        #             current_row = int(x)
        #             if b==len(beamlet_override)-1 and row_counter<CP.num_rows-1:
        #                 break_bool=True
        #                 break
        #
        #             self.left_leaf_position[row_counter] = int(CP.field_position[beamlet_override[b]][1])
        #             print current_row
        #             print self.left_leaf_position
        #             print self.right_leaf_position
        #
        #     # todo fix this conversion to something normal
        #     if not break_bool:
        #         self.right_leaf_position[row_counter] = int(CP.field_position[beamlet_override[-1]][1])+1


        else:
            self.left_leaf_position = aper_left_pos[:]
            self.right_leaf_position = aper_right_pos[:]
            self.intensity = aper_intensity

        self.build_Dkj(CP, data)

    def build_leaf_metadata(self, field, CP):
        self.left_leaf_position = []
        self.right_leaf_position = []

        for row in CP.row_array:

            if len(np.where(field[row][:] > 0)[0]) > 0:
                self.left_leaf_position.append(int(np.where(field[row][:] > 0)[0][0]))
                self.right_leaf_position.append(int(np.where(field[row][:] > 0)[0][-1]))
            else:
                self.left_leaf_position.append(int(field.shape[1] / 2))
                self.right_leaf_position.append(int(field.shape[1] / 2))

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
