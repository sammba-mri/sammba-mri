import os


def _get_dataset_descr(ds_name):
    module_path = os.path.dirname(os.path.abspath(__file__))

    fname = ds_name

    try:
        with open(os.path.join(module_path, 'description', fname + '.rst'))\
                as rst_file:
            descr = rst_file.read()
    except IOError:
        descr = ''

    if descr == '':
        print("Warning: Could not find dataset description.")

    return descr
