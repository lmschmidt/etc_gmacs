# from os.path import dirname,  join

# import numpy as np

from bokeh.plotting import figure
from bokeh.layouts import widgetbox, column, layout, gridplot  # , row
from bokeh.models import ColumnDataSource, glyphs, Div
from bokeh.io import curdoc
from bokeh.models.widgets import Dropdown, RadioButtonGroup, CheckboxButtonGroup, Slider, RangeSlider, Div, Select, Button, TextInput  # Tabs, Panel, Toggle
# from bokeh.embed import components  # for loading stuff,  todo

# import slim as etslim  # reduced python package for gui versions
import strings as stc
import defaults as dfs
import etc as etcalc
# import slim as etslim

# set up data
# create glyphs and their sources
cds_blue = ColumnDataSource(data=dict(xb=[], yb=[]))
cds_red = ColumnDataSource(data=dict(xr=[],  yr=[]))
gly_blue = glyphs.Line(x='xb', y='yb', line_color='blue')
gly_red = glyphs.Line(x='xr', y='yr', line_color='red')

# set up plot

p0 = figure(plot_width=dfs.plot_dims[0],  plot_height=dfs.plot_dims[1], sizing_mode=dfs.plot_sizing_mode,
            x_axis_label=stc.plot_labels[0][1], y_axis_label=stc.plot_labels[0][2])

p0.title.text_font_size='12pt'
p0.xaxis.axis_label_text_font_size = "12pt"
p0.yaxis.axis_label_text_font_size = "12pt"

p0.add_glyph(cds_blue, gly_blue)
p0.add_glyph(cds_red, gly_red)

# !!!!need to impliment https://groups.google.com/a/continuum.io/forum/#!topic/bokeh/nwkfeLbvgUg !!!

# set up widgets
widget_telescope_txt = Div(text="Telescope:")
widget_telescope_size = RadioButtonGroup(labels=stc.telescope_sizes, active=1, name=stc.widget_names[0])
widget_object_type = RadioButtonGroup(labels=stc.object_types, active=0, name=stc.widget_names[1])
widget_star_type = Select(value=stc.star_types_tup[4][0], title=stc.widget_headers[2], options=stc.star_types_tup, name=stc.widget_names[2])
widget_galaxy_type = Select(value=stc.galaxy_types_tup[0][0], title=stc.widget_headers[3], options=stc.galaxy_types_tup, name=stc.widget_names[3])
widget_galaxy_type.disabled = True  # disable on startup
widget_mag_sys = RadioButtonGroup(labels=stc.mag_sys_opts, active=1, name=stc.widget_names[4])
widget_mag = Slider(start=(0), end=(30), value=(25), step=(0.1), title=stc.widget_headers[5], name=stc.widget_names[5],  callback_policy="mouseup")
widget_filter = Select(value=stc.filters_tup[7][0], title=stc.widget_headers[6], options=stc.filters_tup, name=stc.widget_names[6])
widget_grating = RadioButtonGroup(labels=stc.grating_opts, active=0, name=stc.widget_names[7])
widget_moon_txt = Div(text="Days from new moon:")
widget_moon = RadioButtonGroup(labels=stc.moon_opts, active=0, name=stc.widget_names[8])
widget_bin_txt = Div(text="Pixel Binning:")
widget_binning = RadioButtonGroup(labels=stc.bin_opts, active=3, name=stc.widget_names[9])
widget_redshift = Slider(start=(0), end=(10), value=(0), step=(0.01), title=stc.widget_headers[10], name=stc.widget_names[10])
widget_seeing = Slider(start=(0.25), end=(2.0), value=(0.65), step=(0.05), title=stc.widget_headers[11], name=stc.widget_names[11])
widget_slit = Slider(start=(0.25), end=(2.0), value=(0.7), step=(0.05), title=stc.widget_headers[12], name=stc.widget_names[12])
widget_time_inc = RadioButtonGroup(labels=['h','m','s'], active=0, name=stc.widget_names[17])
widget_time = TextInput(value=('1'), title=stc.widget_headers[13], name=stc.widget_names[13])
widget_wavelength = RangeSlider(start=dfs.wavelength_limits[0],  end=dfs.wavelength_limits[1],
                                value=((dfs.wavelength_limits[0]), (dfs.wavelength_limits[1])), step=(10), title=stc.widget_headers[14], name=stc.widget_names[14])
widget_channels = CheckboxButtonGroup(labels=stc.channels,  active=[0, 1],  name=stc.widget_names[15])

widget_plot = Select(value=stc.plot_labels[0][0], title=stc.widget_headers[16], options=[item[0] for item in stc.plot_labels], name=stc.widget_names[16])

widget_update = Button(label="Update", button_type="success")

widget_text = Div(text="Choose input parameters and press the Update button. </br> Controls to Pan, Zoom, Save, or Reset the figure may be selected at the upper right corner of the figure.", width=1000, height=125)

widget_header = Div(text='<h1>'+stc.header1+'</h1><h3>'+stc.header2+'</h3>'+'<hr/>', width=1000, height=100)
widget_footer = Div(text='<hr/>'+'<p>'+stc.footer1+'</p><p>'+stc.footer2+'</p>', width=1000, height=100)

# group widgets for initialization (not layout)
widgets_with_active = [widget_telescope_size, widget_object_type, widget_mag_sys,
                       widget_grating, widget_moon, widget_binning, widget_channels, widget_time_inc]
widgets_with_values = [widget_star_type, widget_galaxy_type, widget_filter, widget_mag,
                       widget_redshift, widget_seeing, widget_slit, widget_time, widget_wavelength, widget_plot]
all_widgets = widgets_with_active + widgets_with_values

# dict to hold all etc input values
etc_inputs = {key: None for key in stc.widget_names}

# link callbacks

def update_etc_inputs(attr, old, new):
    # update dictionary of ETC inputs
    print('updating ETC dict\n')
    # print(etc_inputs, '\n')
    for i, widge in enumerate(widgets_with_values):
        # print(i, widge.name)
        etc_inputs[widge.name] = widge.value

    for i, widge in enumerate(widgets_with_active):
        # print(i, widge.name)
        etc_inputs[widge.name] = widge.active

    # disable object selectors depending on star or galaxy button 
    if etc_inputs['widget_object_type'] == 0:
        widget_galaxy_type.disabled = True
        widget_star_type.disabled = False
    else:
        widget_galaxy_type.disabled = False
        widget_star_type.disabled = True

    # set time label and ranges
    if etc_inputs['widget_time_inc'] == 0:
        widget_time.title = "Exposure Time [h]"
        # widget_time.start = 1
        # widget_time.end = 100
        # widget_time.step = 0.05
    elif etc_inputs['widget_time_inc'] == 1:
        widget_time.title = "Exposure Time [m]"
        # widget_time.start = 1
        # widget_time.end = 60
        # widget_time.step = 1
    elif etc_inputs['widget_time_inc'] == 2:
        widget_time.title = "Exposure Time [s]"
        # widget_time.start = 0.1
        # widget_time.end = 60
        # widget_time.step = 0.1

    # force redshift of 3 for LBG due to template limitations
    if widget_galaxy_type.value == 'lbg_all_flam':
        widget_redshift.value = 3
    
    # ensure that at least one channel is selected
    if not etc_inputs['widget_channels']:
        widget_channels.active = [0,1]

    # Set wavelength range if channels are disabled

    print(etc_inputs)
    
def update_figure():
    # to populate dictionary on first run if no values are changed, 
    # toggle triggers dictionary update, kind of a hack
    if None in etc_inputs.values():
        widget_telescope_size.active = False
        widget_telescope_size.active = True
        #sess = etslim.session(etc_inputs) # create an etc session object with initial values

    print('updating figure \n')
    wavelength, plot_y_red, plot_y_blue, labels, title, messages = etcalc.recalculate(etc_inputs)
    print(labels, title, messages)
    widget_text.update(text='Plotting: ' + title + messages)
    
    # plot_x, plot_yb, plot_yr = sess.update(caller)
    if (etc_inputs['widget_channels'] == [0]):
        cds_blue.data = dict(xb=wavelength, yb=plot_y_blue)
        gly_blue.line_alpha = 0.5
        gly_red.line_alpha = 0.
        p0.title.text = title
        p0.xaxis.axis_label = labels[0]
        p0.yaxis.axis_label = labels[1]
    elif (etc_inputs['widget_channels'] == [1]):
        cds_red.data = dict(xr=wavelength, yr=plot_y_red)
        gly_blue.line_alpha = 0.
        gly_red.line_alpha = 0.5
        p0.title.text = title
        p0.xaxis.axis_label = labels[0]
        p0.yaxis.axis_label = labels[1]
    else:  # crashless catch-all
        gly_blue.line_alpha = 0.5
        gly_red.line_alpha = 0.5
        etc_inputs['widget_channels'] = [0, 1]
        cds_blue.data = dict(xb=wavelength, yb=plot_y_blue)
        cds_red.data = dict(xr=wavelength, yr=plot_y_red)
        p0.title.text = title
        p0.xaxis.axis_label = labels[0]
        p0.yaxis.axis_label = labels[1]

# get current values
for i, widge in enumerate(widgets_with_values):
    # print(i, widge.name)
    # etc_inputs.widge = widge.value
    widge.on_change("value", update_etc_inputs)


for i, widge in enumerate(widgets_with_active):
    # etc_inputs.widge = widge.value
    widge.on_change("active", update_etc_inputs)

widget_update.on_click(update_figure)

# Set up layouts and add to document

widget_group_one = widgetbox(children=[widget_telescope_txt, widget_telescope_size, widget_object_type, widget_star_type, widget_galaxy_type])
widget_group_two = widgetbox(children=[widget_filter, widget_mag, widget_mag_sys])
widget_group_three = widgetbox(children=[widget_grating, widget_redshift, widget_time, widget_time_inc, widget_seeing, widget_slit, widget_moon_txt, widget_moon, widget_wavelength,
                               widget_bin_txt, widget_binning, widget_channels, widget_plot, widget_update], sizing_mode=dfs.plot_sizing_mode)
widget_group_four = column(children=[widget_header, p0, widget_text, widget_footer], sizing_mode=dfs.plot_sizing_mode)
widgets = column(children=[widget_group_one, widget_group_two, widget_group_three], width=dfs.toolbar_width)

curdoc().add_root(layout([[widgets, widget_group_four]],  sizing_mode=dfs.plot_sizing_mode))
curdoc().title = "GMACS ETC 2.0"
