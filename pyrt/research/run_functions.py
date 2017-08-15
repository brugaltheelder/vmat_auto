from pyrt.optimization.vmat import *
from pyrt.tools import print_structure_info


# Run All
def run_all(case,input_dict):
    print_in_box('Running All Cases')
    for filename in os.listdir(input_dict['cwd']):

        if '10' in filename:
            continue

        if '201' in filename:
            continue

        if '202' in filename:
            continue
        if '203' in filename:
            continue

        if '204' in filename:
            continue

        if '205' in filename:
            continue

        if '304' in filename:
            continue

        if case in filename:
            # build input dictionary for a particular case
            input_dict['filename'] = filename

            # call run_case(input)

            run_case(input_dict)
            # print_structure_info(data)
            print '-' * 40
    print_in_box('All Cases Completed')


def run_case(input_dict):

    print_in_box('Now Running {}'.format(input_dict['filename']))
    # build model object
    model = vmat_mip(input_dict)

    # do stuff (optimize)
    model.optimize(run_tag=input_dict['filename'])

    # print out DVH
    model.plot_DVH(run_tag=input_dict['filename'], saveDVH=True, num_bins=500)

    print_in_box('Finished Running {}'.format(input_dict['filename']),1)