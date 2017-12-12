import numpy as np
import scipy.sparse as sps
import re
import matplotlib.pyplot as plt
import seaborn as sns
import itertools



class patient_data(object):
    def __init__(self, input_dict):
        self.input_dict = input_dict.copy()
        self.data_file = self.input_dict['data_file']
        self.beamlet_tags = Read_data_header(path_to_beam_directory+"bixelmap")
        self.beamlet_map = Read_data_img(path_to_beam_directory+"bixelmap",beamlet_tags)

        # these are the voxel assignments to different structures (each structure found in roimask.header, which you can manually type out
        self.roi_index_dict = self.input_dict['roi_index_dict']

        # Mask is 3D array that represents voxel in position (grid with patient somewhere in middle); at each coordimate a value will tell structue assignment for that voxel
        print 'read in mask tags'
        self.mask_tags = Read_data_header(path_to_case_directory + "roimask")
        print 'read in mask img'
        self.mask_map = Read_data_img(path_to_case_directory + "roimask", mask_tags)

        self.priority_voxel_list = self.input_dict['priority_voxel_list']

        for index in priority_voxel_list:
            self.mask_assignment_unique[np.where(flat_mask & 2 ** (index - 1))] = 2 ** (index - 1)
            self.mask_assignment_unique_by_index[np.where(flat_mask & 2 ** (index - 1))] = index

        self.dij_vector, self.bix_vector, self.vox_vector, self.nbix, self.nvox = read_in_dij(path_to_dij)

        # Creat sparse coordinate matrix
        self.dij_coo = sps.coo_matrix((self.dij_vector, (self.bix_vector, self.vox_vector)), shape=(self.nbix, self.nvox))

        # convert to row- or column- sorted
        # self.dij_csr = self.dij_coo.tocsr()
        self.dij_csc = self.dij_coo.tocsc()

        # read in num beamlets and cumulative beamlet thing
        self.beamlet_x, self.beamlet_y, self.num_control_points = beamlet_map.shape

        self.beamlets_per_cp = self.beamlet_x * self.beamlet_y
        self.cumulative_beamlets_per_cp = np.array([0] + np.cumsum(self.beamlets_per_cp).tolist())
        self.num_beamlets = int(self.beamlets_per_cp.sum())

        self.cp_redundancy = self.input_dict['model_params']['cp_redundancy']



        print 'Building Structures'
        self.structures = []
        self.build_structures()

        print 'Building CP'
        self.control_points = []
        self.generate_control_point_data(modality, cp_redundancy=self.cp_redundancy)


    def Read_data_header(file_tag, interested_tags=['x_dim', 'y_dim', 'z_dim', 'data_type', 'byte_order']):
        header_filename = file_tag + ".header"
        img_filename = file_tag + ".img"
        return_dict = {}
        with open(header_filename, 'r')as f:
            for line in f.readlines():
                a = re.sub(r'[\n\r]+', '', line).split(" ")
                if len(a) == 3:
                    if a[0] in interested_tags:
                        return_dict[a[0]] = a[2]
        return return_dict

    def Read_data_img(file_tag, input_dict):
        img_filename = file_tag + ".img"
        dt = np.float32
        if input_dict['data_type'] == 'float':
            dt = np.float32
        elif input_dict['data_type'] == 'uint32':
            dt = np.uint32
        else:
            print "data type input issues"
        order = 'C'
        if input_dict['byte_order'] == '0':
            order = 'F'
        # Read in 1D image array (block of memory) and convert bits into float value and put it back into 3D
        input_img = np.fromfile(img_filename, dtype=dt)
        # print input_img.shape
        input_img = input_img.reshape((int(input_dict['x_dim']), int(input_dict['y_dim']), int(input_dict['z_dim'])),
                                      order=order)
        # print input_img.shape

        return input_img

    # here is an example function to get organ assignment per voxel (you can make this faster/with arrays...but here is an example)
    def get_voxel_assignment(vox_mask_value, vox_structure_dict):
        structures_belonging = []
        for index, structure in vox_structure_dict.iteritems():
            if vox_mask_value & (2 ** (index - 1)):
                structures_belonging.append(structure)
        if len(structures_belonging) == 0:
            return ['Air']
        else:
            return structures_belonging

    # here is an example function to get organ assignment per voxel (you can make this faster/with arrays...but here is an example)
    def belongs_to_structure(vox_mask_value, structure_index):
        if vox_mask_value & (2 ** (structure_index - 1)):
            return True
        return False

    def read_in_dij(filepath):

        with open(filepath + 'Size_out.txt', 'r')as f:
            lines = f.readlines()
            nvox = int(lines[0].strip())
            nbix = int(lines[1].strip())
            ndij = int(lines[2].strip())

        bix_vector = np.fromfile(open(filepath + 'Bixels_out.bin', 'r'), dtype=np.uint32)
        vox_vector = np.fromfile(open(filepath + 'Voxels_out.bin', 'r'), dtype=np.uint32)
        dij_vector = np.fromfile(open(filepath + 'Dijs_out.bin', 'r'), dtype=np.float32)

        # print vox_vector.shape
        return (dij_vector, bix_vector, vox_vector, nbix, nvox)


    def build_structures(self):
        # Create list of real structure names
        structure_names = []
        structure_sizes = {}
        structure_index = {}

        for index, structure in self.roi_index_dict.iteritems():
            print index, structure
            structure_names.append(structure)
            structure_index[structure_names[-1]] = index
            print 'finding structure voxel indices for structure: {}'.format(structure)
            structure_sizes[structure_names[-1]] = np.where(self.mask_assignment_unique_by_index == index)[0]
            # if num_voxels[structure_names[-1]].size > 0:
            #     print '{} voxels in structure {}'.format(num_voxels[structure_names[-1]].size, structure)


        # gather init data for struct
        # Get all volume/structure/data names

        flattened_beamlet_map = beamlet_map.flatten(order='F')

        for s in range(structure_names.size):

            if structure_sizes[s]>0:
                # set prescription, is_target
                Rx, is_target = 0., False

                if name in self.input_dict['Rx'].keys():
                    Rx = self.input_dict['Rx'][name]
                    is_target = True

                print 'slicing Dijs for structure {}'.format(s.name)
                Dijs = self.dij_csc[:, structure_sizes[s]]
                Dijs_active_per_structure = []


                starting_index = 0
                ending_index = (int(self.num_control_points) + 1) * beamlets_per_beam  # non inclusive
                indices = starting_index + np.where(flattened_beamlet_map[starting_index:ending_index] == 1)[0]

                Dijs_active_per_structure.append(Dijs[indices, :])





                self.structures.append(structure(name=name,index=structure_index[name], f=self.f, Rx=Rx, num_vox=structure_sizes[name],
                                                 num_beamlets=self.num_beamlets, is_target=is_target, data_file=self.data_file))










































