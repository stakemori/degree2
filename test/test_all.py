import os
import re
import imp

from os.path import dirname, isfile, join
pat = re.compile(".*test.+py$")
this_dir = dirname(__file__)

def test_module_names():
    return [f for f in os.listdir(this_dir) if re.match(pat, f)
            and isfile(join(this_dir, f))]

def import_tests():
    for f in test_module_names():
        mod_name = "degree2.test." + f.split(".")[0]
        pth = join(this_dir, f)
        imp.load_source(mod_name, pth)

import_tests()
