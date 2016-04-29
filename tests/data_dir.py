import os.path as opath
from sage.all import load


def data_dir(name):
    return opath.join(opath.dirname(opath.abspath(__file__)), "data", name)


def load_from_data_dir(fname, dirname):
    return load(opath.join(data_dir(dirname), fname))
