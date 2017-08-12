__author__ = 'troy'


def print_in_box(string_input):
    length = len(string_input)
    print '-'*(length+6)
    print '|  {}  |'.format(string_input)
    print '-'*(length+6)