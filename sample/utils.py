import os
import yaml

def read_yaml(yaml_file):
    with open(yaml_file, 'r') as stream:
        return yaml.safe_load(stream)

def make_dir(d):
    '''
    Makes dir if it does not already exist
    '''

    if not os.path.exists(d):
        os.makedirs(d)

    return

def get_base_dir():
    '''
    base dir is parent repo dir
    '''
    module_dir = get_module_dir()
    return os.path.dirname(module_dir)

def get_module_dir():
    return os.path.dirname(__file__)

def get_plot_dir():
    return os.path.join(get_base_dir(), 'plots')

def get_out_dir():
    return os.path.join(get_base_dir(), 'out')

def get_data_dir():
    return os.path.join(get_base_dir(), 'data')

def get_test_dir():
    module_dir = get_module_dir()
    return os.path.join(module_dir, 'tests')

BASE_DIR = get_base_dir()
MODULE_DIR = get_module_dir()
TEST_DIR = get_test_dir()
PLOT_DIR = get_plot_dir()
OUT_DIR = get_out_dir()
DATA_DIR = get_data_dir()
