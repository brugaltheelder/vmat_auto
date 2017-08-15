import matplotlib.pyplot as plt
import numpy as np

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

def plot_DVH(model, saveName='', showPlots=False, saveDVH=False, run_tag=None,num_bins = 500):
    if run_tag is not None:
        dose_per_struct = [np.copy(model.dose_dict[run_tag][s]) for s in range(len(model.data.structures))]
        run_label = run_tag
    else:
        dose_per_struct = [np.copy(model.current_dose_per_structure[s]) for s in range(len(model.data.structures))]
        run_label='current'

    plt.clf()
    for s in range(len(model.data.structures)):

        hist, bins = np.histogram(dose_per_struct[s], bins=num_bins)
        dvh = 1. - np.cumsum(hist) / float(model.data.structures[s].num_vox)
        plt.plot(bins[:-1], dvh, label=model.data.structures[s].name, linewidth=2)
    lgd = plt.legend(fancybox=True, framealpha=0.5, bbox_to_anchor=(1.05, 1), loc=2)

    plt.title('DVH for run tag: {}'.format(run_tag))

    plt.xlabel('Dose')
    plt.ylabel('Fractional Volume')

    plt.gca().set_xlim(left=0.)
    plt.gca().set_ylim(bottom=0., top=1.)


    if len(saveName) > 1 or saveDVH:

        plt.savefig(model.data.input_dict['figure_directory'] + 'dvh_' + saveName + '_' + run_label + '.png',
                    bbox_extra_artists=(lgd,), bbox_inches='tight')
    if showPlots:
        plt.show()

def plot_fluence(model,saveName='', showPlots=False, saveDVH=False, run_tag=None, aperture_list = []):
    pass

def plot_fluence_map(fluence):

    pass

def is_adjacent(aper_L, aper_R, distance_limit):
    # for each row, check if L and R leafs are within distance limit
    pass




