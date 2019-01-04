#!/usr/bin/python3
from pathlib import Path

''' paths '''
# get parent directory
parentDir = Path.cwd()  / 'etc_gmacs'

#doc file
doc_path = parentDir / 'readme.pdf'

galaxy_path = parentDir / 'core' / 'kinney'
stellar_path = parentDir / 'core' / 'pickle'
skyfiles_path = parentDir / 'core' / 'skybackground'
filter_path = parentDir / 'core' / 'filters'
atmo_ext_path = parentDir / 'core' / 'atmo_extinction.dat'
dichroic_path = parentDir / 'core' / 'dichroic5577.txt'
grating_path = [parentDir / 'core' / 'grating_red_low_315-1100.txt',  parentDir / 'core' / 'grating_blue_low_315-1100.txt']
ccd_path = [parentDir / 'core' / 'e2v-astro-multi-2-DD.txt',  parentDir / 'core' / 'e2v_blue.txt']
vega_file = parentDir / 'core' / 'alpha_lyr_stis_005.ascii'
mirror_file = parentDir / 'core/mirror_300-1200nm.txt'
galaxy_files = ['SB1', 'SB2', 'SB3', 'SB4', 'SB5', 'SB6', 'S0', 'Sa', 'Sb', 'Sc', 'bulge', 'ellipticals', 'lbg_all_flam']
galaxy_labels = ['Starbursts with E(B-V) < 0.10', 'Starbursts with  0.11 < E(B-V) < 0.21', 'Starbursts with  0.25 < E(B-V) < 0.35',
                 'Starbursts with  0.39 < E(B-V) < 0.50', 'Starbursts with  0.51 < E(B-V) < 0.60', 'Starbursts with  0.61 < E(B-V) < 0.70',
                 'S0 Galaxies', 'Sa Galaxies', 'Sb Galaxies ', 'Sc Galaxies', 'Galactic Bulges', 'Elliptical Galaxies', 'Lyman-break Galaxies']
stellar_files = ['o5v.dat', 'b0v.dat', 'b57v.dat', 'a0v.dat', 'a5v.dat', 'f0v.dat', 'f5v.dat', 'g0v.dat', 'g5v.dat', 'k0v.dat', 'k5v.dat', 'm0v.dat', 'm5v.dat']
skyfiles = ['00d_315-1200nm.csv', '03d_315-1200nm.csv', '07d_315-1200nm.csv', '10d_315-1200nm.csv', '14d_315-1200nm.csv']
filter_files = ['BesU.dat', 'BesB.dat', 'BesV.dat', 'BesR.dat', 'BesI.dat', 'u.dat', 'g.dat', 'r.dat', 'i.dat', 'z.dat']
# filter_files = ['photonUX.dat', 'photonB.dat', 'photonV.dat', 'photonR.dat', 'photonI.dat', 'u.dat', 'g.dat', 'r.dat', 'i.dat', 'z.dat']
