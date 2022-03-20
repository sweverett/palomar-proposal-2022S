import os

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

def get_data_dir():
    return os.path.join(get_base_dir(), 'data')

def get_test_dir():
    module_dir = get_module_dir()
    return os.path.join(module_dir, 'tests')

BASE_DIR = get_base_dir()
MODULE_DIR = get_module_dir()
TEST_DIR = get_test_dir()
