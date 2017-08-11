import numpy as np
import scipy.sparse as sps
import os
import h5py


class patient_data(object):
    def __init__(self, input_dict):
        self.f = h5py.File(input_dict['cwd'] + input_dict['filename'], 'r')
        self.input_dict = input_dict.copy()
        self.structures = []
        self.build_structures()


    def build_structures(self):
        # Create list of real structure names
        self.structure_names = []
        structure_sizes = []
        patient_structure_names = self.f['patient/StructureNames']
        patient_structure_sizes = self.f['patient/SampledVoxels']
        for i in range(patient_structure_names.size):
            self.structure_names.append(''.join(chr(j) for j in self.f[patient_structure_names[i][0]][:]))
            structure_sizes.append(int(self.f[patient_structure_names[i][0]].shape[0]))

        #gather init data for struct
        #Get all volume/structure/data names
        data_matrix = self.f['data/matrix']
        for s in range(data_matrix['A'].size):

            A_ref = data_matrix['A'][s]

            self.structures.append(structure())
            self.structure_names.append(structure(self.structure_names[i], self.f, self.input_dict, i))



class structure(object):
    def __init__(self, name, f, input_dict, i):
        self.name = name

        # Access Matrix subgroup in patient data (Dose to points, structure name, etc.)
        # f['data/matrix'].keys()

        # Determine corresponding number of voxels for structure
        c = f['patient/SampledVoxels']
        self.structure_size = f[c[i][0]].shape[0]
        # print name,self.structure_size






        # patient_structure_size = {}
        # for i in range(c.size):
        #     patient_structure_size[structure_names_real[i]] = f[c[i][0]].shape[0]
        # # print patient_structure_size


        # # get ref to data/matrix
        # b = f['data/matrix']
        #
        # # read in names
        # structure_names = []
        # for i in range(b['Name'].size):
        #     structure_names.append(''.join(chr(j) for j in f[b['Name'][i][0]][:]))
        # # print structure_names
        # print '-' * 40
        #
        # # read in A's
        # A_mat = {}
        # for i in range(b['A'].size):
        #     if structure_names[i] not in structure_names_real:
        #         continue
        #     else:
        #         print i,
        #         # Find sparse matrices/structures:
        #         if np.asarray(f[b['A'][i][0]]).shape == (3,):
        #             print 'SPARSE'
        #             print structure_names[i]
        #             indices = np.asarray(f[b['A'][i][0]]['jc'])
        #             indptr = np.asarray(f[b['A'][i][0]]['ir'])
        #             data_values = np.asarray(f[b['A'][i][0]]['data'])
        #             print 'now testing structure {}'.format(structure_names[i])
        #             print 'this should have a shape of {} by {}'.format(patient_structure_size[structure_names[i]],
        #                                                                 indices.shape[0] - 1)
        #             # csr_matrix((data, indices, indptr), [shape=(M, N)])
        #             A_mat[structure_names[i]] = sps.csr_matrix((data_values, indptr, indices), shape=(
        #             indices.shape[0] - 1, patient_structure_size[structure_names[i]]))
        #             #           print A_mat[structure_names[i]].size, A_mat[structure_names[i]].shape
        #             np.setdiff1d(np.arange(1, 11788), np.unique(indptr)).shape
        #
        #         else:
        #             A_mat[structure_names[i]] = sps.csr_matrix(np.asarray(f[b['A'][i][0]]))
        #
        # print ''
        # print '-' * 40

#
# #Access Misc subgroup in patient data
# f['data/misc'].keys()
#
# # get ref to data/matrix
# c = f['data/misc']
#
# #Assign Rx to structure
# Rx={}
# for i in range(b['A'].size):
#     if structure_names[i] not in structure_names_real:
#         continue
#     elif str(structure_names[i])=='PTV':
#         Rx[structure_names[i]]= c['InitialiseReferenceDose'][0][0]
#     elif str(structure_names[i])=='PTV Vesicles':
#         Rx[structure_names[i]]= c['InitialiseReferenceDose'][1][0]
#     else:
#         Rx[structure_names[i]]=0
#
# print Rx
#
# #Access Misc subgroup in patient data
# f['data/misc'].keys()
#
# # get ref to data/matrix
# b = f['data/misc']
#
# print b['InitialiseReferenceDose'].size
# print b['InitialiseReferenceDose'][1][0]
#
# print ''
# print '-'*40
#
#
# # cwd = '/Users/jtmargo/Desktop/Josh/Clemson University/Research/Treatment Plan/TROTS Data/Prostate_VMAT/'
# # for filename in os.listdir(cwd):
# #     print filename
# #     f = h5py.File(cwd+filename,'r')
# #     #     get ref to data/matrix
# #     b = f['data/misc']
# #     print b['InitialiseReferenceDose'].size, np.asarray(b['InitialiseReferenceDose'])
#
#
