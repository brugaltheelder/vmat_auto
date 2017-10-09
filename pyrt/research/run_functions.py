from pyrt.optimization.vmat import *
from pyrt.tools import print_structure_info


# Run All
def run_all(case,input_dict):
    print_in_box('Running All Cases')

    if case == 'Prostate_VMAT':
        skip_files = [210,304,305]
    else:
        skip_files = []

    if case == "Head-and-Neck" :
        skip_files = [01,02,04,06,07]
    else:
        skip_files=[]

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

    print_in_box('Now Running {}'.format(input_dict['filename'][0:-4]),1)
    # build model object
    model = vmat_mip(input_dict)

    # do stuff (optimize)
    model.optimize(run_tag=input_dict['filename'][0:-4])

    # print out DVH
    model.plot_DVH(run_tag=input_dict['filename'][0:-4], saveDVH=True, num_bins=500, specific_directory =model.data.input_dict['case_directory'])

    #Print out aperture shapes
    plot_all_selected_apertures(model)

    print_in_box('Finished Running {}'.format(input_dict['filename'][0:-4]),1)