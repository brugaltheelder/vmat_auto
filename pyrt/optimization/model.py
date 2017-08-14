# this file will have "model" structured files

from time import time

import matplotlib.pyplot as plt
import scipy.optimize as spo

from pyrt.data.data_trots import *
from pyrt.tools import *


class model_base(object):
    def __init__(self, input_dict, modality='undefined'):
        self.modality = modality
        print_in_box('Reading in data')
        self.read_in_data(input_dict)
        print_in_box('Data reading completed')
        self.current_dose_per_structure = [np.zeros(s.num_vox) for s in self.data.structures]
        self.current_obj = 0.
        self.dose_dict = {}
        self.obj_dict = {}

    def read_in_data(self,input_dict):
        self.data = patient_data(input_dict)

    def generate_optimization_variables(self):
        print 'Implement in child classes'
        pass

    def optimize(self):
        print 'Implement in child classes'
        pass

    def calc_dose_from_variables(self,x=None):
        print 'Implement in child classes'
        pass

    def save_current_dose(self,run_tag):
        self.dose_dict[run_tag] = [np.copy(self.current_dose_per_structure[s]) for s in range(len(self.data.structures))]


    def plot_DVH(self, saveName='', showPlots=False, saveDVH=False, run_tag=None,num_bins = 100):
        if run_tag is not None:
            dose_per_struct = [np.copy(self.dose_dict[run_tag][s]) for s in range(len(self.data.structures))]
            run_label = run_tag
        else:
            dose_per_struct = [np.copy(self.current_dose_per_structure[s]) for s in range(len(self.data.structures))]
            run_label='current'

        plt.clf()
        for s in range(len(self.data.structures)):

            hist, bins = np.histogram(dose_per_struct[s], bins=num_bins)
            dvh = 1. - np.cumsum(hist) / float(self.data.structures[s].num_vox)
            plt.plot(bins[:-1], dvh, label=self.data.structures[s].name, linewidth=2)
        lgd = plt.legend(fancybox=True, framealpha=0.5, bbox_to_anchor=(1.05, 1), loc=2)

        plt.title('DVH for run tag: {}'.format(run_tag))

        plt.xlabel('Dose')
        plt.ylabel('Fractional Volume')

        plt.gca().set_xlim(left=0.)
        plt.gca().set_ylim(bottom=0., top=1.)


        if len(saveName) > 1 or saveDVH:

            plt.savefig(self.data.input_dict['figure_directory'] + 'dvh_' + saveName + '_' + run_label + '.png',
                        bbox_extra_artists=(lgd,), bbox_inches='tight')
        if showPlots:
            plt.show()












