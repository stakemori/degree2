# -*- coding: utf-8; mode: sage -*-
from sage.all import load, QQ
from os import path
current_dir = path.dirname(path.abspath(__file__))
cached_dir = path.join(current_dir, "cached_data")

from degree2.deg2_fourier import Deg2ModularFormQseries, Deg2EisensteinQseries, KlingenEisenstein_and_cuspforms_space
import degree2.deg2_fourier
_deg2fc_gens_dict = load(path.join(cached_dir, '_fc_dict34.sobj'))

es4 = Deg2ModularFormQseries(4, _deg2fc_gens_dict[4], 34)
degree2.deg2_fourier.Deg2global_gens_dict['es4'] = es4

es6 = Deg2ModularFormQseries(6, _deg2fc_gens_dict[6], 34)
degree2.deg2_fourier.Deg2global_gens_dict['es6'] = es6

x10 = Deg2ModularFormQseries(10, _deg2fc_gens_dict[10], 34)
degree2.deg2_fourier.Deg2global_gens_dict['x10'] = x10

x12 = Deg2ModularFormQseries(12, _deg2fc_gens_dict[12], 34)
degree2.deg2_fourier.Deg2global_gens_dict['x12'] = x12

x35 = Deg2ModularFormQseries(35, _deg2fc_gens_dict[35], 34)
degree2.deg2_fourier.Deg2global_gens_dict['x35'] = x35

load(path.join(cached_dir, "tuples.py"))
