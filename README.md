GMACS Exposure Time Calculator

etc_gmacs is an update to the original GMACS exposure time calculator.  The original calculator was a combination of a .cgi script web form with the calculations performed by Octave, an open source Matlab equivalent. Several of the packages were very old and several desired features had not yet been implimented. 

This version is written for Python 3.7 and uses [Bokeh](https://bokeh.pydata.org/en/latest/) to generate the etc input widgets and interactive plot.

Packages required (use pip install if not already installed):

bokeh

astropy

pandas

pathlib

spectres

matplotlib

scipy

To RUN:

1.  After cloning the repository you should have a folder called etc_gmacs.  
2.  Open up a command prompt (this was developed in Windows using Anaconda, so open an Anaconda prompt, will likely work in linux as well, but have not tested).
3.  Navigate to one level above etc_gmacs (for example if the project folder is located at /Documents/etc_gmacs, run the following command from the /Documents folder)
4.  run bokeh serve etc_gmacs
5.  You should see some messages about starting the server, then you can use your web browser to navigate to the given address (probably something like http://localhost:5006/etc_gmacs)
6.  ETC should be ready to use!

In typical fashion, this was coded quickly, so there are possibly things that don't quite work right, code that is hard to follow, etc. If you notice any errors, please contact me.
