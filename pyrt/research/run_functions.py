from pyrt.optimization.vmat import *
from pyrt.tools import print_structure_info


# Run All
def run_all(case,input_dict):
    print_in_box('Running All Cases')

    if case == 'Prostate_VMAT':
        skip_files = [304,305]
    else:
        skip_files = []

    for filename in os.listdir(input_dict['cwd']):




        if case in filename:
            if int(filename[-7:-4]) not in skip_files:
                continue
            # build input dictionary for a particular case
            input_dict['filename'] = filename

            # call run_case(input)

            run_case(input_dict)
            print '-' * 40
    print_in_box('All Cases Completed')


def run_case(input_dict):

    print_in_box('Now Running {}'.format(input_dict['filename']),1)
    # build model object
    model = vmat_mip(input_dict)

    # do stuff (optimize)
    model.optimize(run_tag=input_dict['filename'])

    # print out DVH
    model.plot_DVH(run_tag=input_dict['filename'], saveDVH=True, num_bins=500)

    #Print out aperture shapes
    plot_all_selected_apertures()

    print_in_box('Finished Running {}'.format(input_dict['filename']),1)