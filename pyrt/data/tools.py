
import numpy as np
import scipy.sparse as sps
import h5py




def find_min_max_row(data):
    pass



class structure(object):
    def __init__(self, name,index,  A_ref, f, Rx, num_vox, num_beamlets, is_target=False):
        self.name = name
        self.index = index
        self.rx = Rx
        self.num_vox = num_vox
        self.Dij = None
        self.num_beamlets = num_beamlets
        self.is_target = is_target
        self.import_dose(A_ref, f)

    def import_dose(self, A_ref, f):
        if np.asarray(f[A_ref[0]]).shape == (3,):
            print 'importing {} Dij as sparse matrix'.format(self.name)
            indices = np.asarray(f[A_ref[0]]['jc'])
            indptr = np.asarray(f[A_ref[0]]['ir'])
            data = np.asarray(f[A_ref[0]]['data'])
            # sanity check
            if self.num_beamlets != indices.size - 1:
                print 'ERROR IN DIMENSION MISMATCH' * 40

            self.Dij = sps.csr_matrix((data, indptr, indices), shape=(self.num_beamlets, self.num_vox))
        else:
            print 'importing {} Dij as dense matrix, converting to sparse...'.format(self.name)
            self.Dij = sps.csr_matrix(np.asarray(f[A_ref[0]]))


class control_point(object):

    def __init__(self, cp_number, field, min_row, max_row, initial_beamlet_index, number_beamlets ):
        self.cp_number = cp_number
        self.min_row = min_row
        self.max_row = max_row
        self.number_beamlets = number_beamlets
        self.initial_beamlet_index = initial_beamlet_index
        self.final_beamlet_index = initial_beamlet_index+number_beamlets

        self.build_leaf_metadata(field)
        # check if beamlet row sum == number_beamlets

    def build_leaf_metadata(self,field):
        pass