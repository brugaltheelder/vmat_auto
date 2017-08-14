
import numpy as np

class aperture(object):
    def __init__(self, CP, aperture_shape_details):
        self.left_leaves = None
        self.right_leaves = None
        self.intensity = None
        self.Dkj_per_structure = None


    def calc_dose(self,x=None):
        pass