# Run All
def run_all(input_dict):
    for filename in os.listdir(cwd):
        # build input dictionary for a particular case
        input_dict['filename'] = filename
        # call run_case(input)
        run_case(input_dict)
        print_structure_info(data)
        print '-' * 40


def run_case(input_dict):

    # build model object
    model = vmat_mip(input_dict)

    # do stuff (optimize)
    model.optimize()

    # print out DVH
    model.plot_DVH(run_tag='default', saveDVH=True, num_bins=500)