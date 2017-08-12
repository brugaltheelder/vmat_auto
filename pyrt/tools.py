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