
import numpy as np
import scipy.sparse as sps
import h5py




def find_min_max_row(data):
    b = data.f['patient/Beams/BeamConfig']
    min_row = np.asarray(data.f[b['Field'][0][0]]).shape[1] # [np.asarray(data.f[b['Field'][0][0]]).shape[1]]*data.num_control_points
    max_row = 0 # [0]*data.num_control_points


    for cp in range(data.num_control_points):
        field = np.asarray(data.f[b['Field'][cp][0]])

        min_cp = int(np.where(field == 1)[0])
        max_cp = int(np.where(field == np.asarray(data.f['patient/Beams/ElementIndex'][0][cp]))[0])
        if min_cp < min_row:
            min_row = min_cp
        if max_cp > max_row:
            max_row = max_cp

    return min_row, max_row+1


def find_min_max_row_imrt(data):
    b = data.f['patient/Beams/BeamConfig']
    min_row = [np.asarray(data.f[b['Field'][0][0]]).shape[1]]*data.num_control_points
    max_row = [0]*data.num_control_points


    for cp in range(data.num_control_points):
        field = np.asarray(data.f[b['Field'][cp][0]])

        min_row[cp] = int(np.where(field == 1)[0])
        max_row[cp]= int(np.where(field == np.asarray(data.f['patient/Beams/ElementIndex'][0][cp]))[0])+1


    return min_row[:], max_row[:]


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


class control_point_vmat(object):

    def __init__(self, cp_number, field, min_row, max_row, initial_beamlet_index, number_beamlets ):
        self.cp_number = cp_number
        self.min_row = min_row
        self.max_row = max_row
        self.num_rows = max_row-min_row
        self.number_beamlets = number_beamlets
        self.initial_beamlet_index = initial_beamlet_index
        self.final_beamlet_index = initial_beamlet_index+number_beamlets
        self.field = field.copy()

        self.build_leaf_metadata(field)
        self.build_spatial_metadata(field)

        # check if beamlet row sum == number_beamlets
        if np.array(self.width_per_row).sum() != self.number_beamlets:
            print np.array(self.width_per_row).sum(), self.number_beamlets
            print "ERROR: NOT ALL BEAMLETS COUNTED or GAPS IN FLUENCE MAP"

    def build_spatial_metadata(self,field):
        self.field_position =  np.zeros((self.number_beamlets,2))
        beamlet_counter = 0
        for r in range(self.num_rows):
            for i in range(self.width_per_row[r]):
                self.field_position[beamlet_counter,0] = int(r+self.min_row)
                self.field_position[beamlet_counter,1] = int(i+self.left_leaf_position[r])

                beamlet_counter+=1




    def build_leaf_metadata(self,field):
        self.left_leaf_position = []
        self.left_leaf_index = []
        self.width_per_row = []
        self.row_array = range(self.min_row, self.max_row)

        for row in self.row_array:

            self.left_leaf_position.append(int(np.where(field[row][:] > 0)[0][0]))
            self.left_leaf_index.append(field[row][np.where(field[row][:] > 0)[0][0]] - 1)
            self.width_per_row.append(int(np.argmax(field[row][:])) - int(np.where(field[row][:] > 0)[0][0]) + 1)
            # else:
            #     self.left_leaf_position.append(0)
            #     self.left_leaf_index.append(-1)
            #     self.width_per_row.append(0)





