import os
import re
# import imp

from os.path import dirname, isfile, join
pat = re.compile(".*test.+py$")
this_dir = dirname(__file__)


def test_module_names():
    return [f for f in os.listdir(this_dir) if re.match(pat, f)
            and isfile(join(this_dir, f)) and join(this_dir, f) != __file__]

from degree2.tests import (test_fc_mul_add, test_eigenforms, test_divide,
                           test_save_load, test_interpolate,
                           test_gens, misc_test, test_prec_class,
                           test_vector_valued, test_const)

# def import_tests():
#     for f in test_module_names():
#         mod_name = "degree2.tests." + f.split(".")[0]
#         pth = join(this_dir, f)
#         imp.load_source(mod_name, pth)

# import_tests()
