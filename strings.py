#!/usr/bin/python3

import paths as etpaths

# title,  x,  y
plot_labels = [('Signal-to-Noise Ratio', 'Wavelength [\u212b]', 'SNR px\u207b\u00b9'),
               ('Observed Spectrum', 'Wavelength [\u212b]', 'Counts px\u207b\u00b9'),
               ('Observed Spectrum + Noise', 'Wavelength [\u212b]', 'Counts px\u207b\u00b9'),
               ('Observed Sky Background', 'Wavelength [\u212b]', 'Counts px\u207b\u00b9'),
               ('Dichroic Throughput', 'Wavelength [\u212b]', 'Throughput'),
               ('Grating Efficiency', 'Wavelength [\u212b]', 'Efficiency'),
               ('CCD Quantum Efficiency', 'Wavelength [\u212b]', 'Quantum Efficiency'),
               ('Atmospheric Extinction', 'Wavelength [\u212b]', 'Extinction')]

widget_headers = ["Telescope Mode",  # 0
                  "Object Type",  # 1
                  "Stellar Classification",  # 2
                  "Galactic Classifications",  # 3
                  "Magnitude System",  # 4
                  "Magnitude",  # 5
                  "Filter",  # 6
                  "Grating",  # 7
                  "Days since/until new moon",  # 8
                  "Pixel Binning",  # 9
                  "Redshift [z]",  # 10
                  "Seeing [arcsec]",  # 11
                  "Slit Width [arcsec]",  # 12
                  "Exposure Time [h]",  # 13
                  "Spectral Range [\u212b]",  # 14
                  "Active Channels",  # 15
                  "Plot Type"]  # 16

widget_names = ['widget_telescope',
                'widget_object_type',
                'widget_star_type',
                'widget_galaxy_type',
                'widget_mag_sys',
                'widget_mag',
                'widget_filter',
                'widget_grating',
                'widget_moon',
                'widget_binning',
                'widget_redshift',
                'widget_seeing',
                'widget_slit',
                'widget_time',
                'widget_wavelength',
                'widget_channels',
                'widget_plot',
                'widget_time_inc']  # matched w/ header nums except for time_inc as no header is needed

telescope_sizes = ["First light", "Full Size"]
object_types = ["Stellar", "Galactic"]

header1 = 'GMACS : Exposure Time Calculator v2.0'
header2 = 'Munnerlyn Astronomical Instrumentation Lab'
footer1 = '<a href='+str(etpaths.doc_path)+'>README</a>'
footer2 = '''Send questions or bug reports to <a href='&#109;&#97;i&#108;to&#58;%6Cs&#99;h%6Di&#100;t&#64;ph%&#55;9&#115;%&#54;9&#99;&#115;&#46;%74amu&#46;%65du'>lsch&#109;idt&#64;p&#104;y&#115;ics&#46;t&#97;m&#117;&#46;&#101;du</a>.'''
mag_sys_opts = ['Vega', 'AB']
grating_opts = ['Low Resolution', 'High Resolution']
filter_opts = [name[:-4] for name in etpaths.filter_files]
moon_opts = ['0', '3', '7', '10', '14']
bin_opts = ['1x1', '2x2', '3x3', '4x4', 'Nyquist', 'Resel']
noise_opts = ['With Noise', 'Without Noise']
channels = ['Blue Channel', 'Red Channel']

star_types = [name[:-4] for name in etpaths.stellar_files]
galaxy_types = etpaths.galaxy_files  # these files lack suffixes
galaxy_labels = etpaths.galaxy_labels
star_types_tup = [(i, i.upper()) for i in star_types]  # unless there's a need to distinguish...
galaxy_types_tup = list(zip(galaxy_types,galaxy_labels))
filters_tup = [(etpaths.filter_files[i], filter_opts[i]) for i in range(len(filter_opts))]
