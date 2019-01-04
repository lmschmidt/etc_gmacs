GMACS Exposure Time Calculator

etc_gmacs is an update to the original GMACS exposure time calculator.  The original calculator was a combination of a .cgi script web form with the calculations performed by Octave, an open source Matlab equivalent. Several of the packages were very old and several desired features had not yet been implimented. 

This version is written for Python 3.7 and uses [Bokeh](https://bokeh.pydata.org/en/latest/) to generate the etc input widgets and interactive plot.

## Start Bokeh Server
The Bokeh server runs on the webserver and is visible from the localhost.  A reverse proxy forwards external requests to the local host to make the server visible to the outside world.  To start the ETC, in a ssh session run the following command:
```
nohup bokeh serve etc_gmacs --port 5100 --allow-websocket-origin=instrumentation.tamu.edu &
```
https://superuser.com/questions/448445/run-bash-script-in-background-and-exit-terminal

## Stop Bokeh Server
To stop the server in order to update, ssh to the webserver and traverse to the etc_gmacs folder. Run the following
```
ps -e | grep bokeh
```
This will list the PID of the server process which you can then kill, 
```
kill xxxxxx
```
Where xxxxxx is the PID.
