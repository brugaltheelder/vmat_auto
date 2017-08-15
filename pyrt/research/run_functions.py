from pyrt.optimization.vmat import *
from pyrt.tools import print_structure_info


# Run All
def run_all(case,input_dict):
    for filename in os.listdir(cwd):

        if case in filename:
            # build input dictionary for a particular case
            input_dict['filename'] = filename

            # call run_case(input)
            run_case(input_dict)
            # print_structure_info(data)
            print '-' * 40


def run_case(input_dict):

    # build model object
    model = vmat_mip(input_dict)

    # do stuff (optimize)
    model.optimize()

    # print out DVH
    model.plot_DVH(run_tag=input_dict['filename'], saveDVH=True, num_bins=500)