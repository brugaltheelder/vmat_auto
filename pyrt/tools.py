__author__ = 'troy'


def print_in_box(string_input):
    length = len(string_input)
    print '-'*(length+6)
    print '|  {}  |'.format(string_input)
    print '-'*(length+6)



def print_structure_info(data):
    print '-'*20
    for s in data.structures:
        print s.name,s.rx, s.num_vox, s.num_beamlets, s.is_target, s.Dij.shape
        print '-'*20



def plot_fluence_map(fluence,field):
    pass

def is_adjacent(aper_L, aper_R, distance_limit):
    # for each row, check if L and R leafs are within distance limit
    pass




class aperture(object):
    def __init__(self, CP, aperture_shape_details):
        self.left_leaves = None
        self.right_leaves = None
        self.intensity = None
        self.Dkj_per_structure = None


    def calc_dose(self,x=None):
        pass