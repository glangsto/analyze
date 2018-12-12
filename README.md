## analyze
### analyize: plotting tools for reviewing/averaging/calibrated spectra line observations from Science Aficionado Telescopes.   

The plotting programs expect Ascii input spectra with 
extensive headers describing the observations.

The plotting programs are compatible with spectra created using any of the GnuRadio NsfIntegrate??.grc designs
These python programs were created using the Gnu Radio Companion (GRC).   These programs
contained here are for analysis after the observations, not for data taking.

These programs are also compatible with the older NsfWatch output (watch.py)

### Programs:

* S     Summarize a set of observations.  Function calls s.py
* C     _Main_ program which calibrates, averages, flags and plots sets of spectra
* R     Function to call r.py, which plots raw spectra

### Python:

* c.py  Large python function to calibrate, average, flag and baseline fit observations.
* r.py  Python function to read and plot the _raw_ spectra
* s.py  Python function to read all selected spectra in a directory and summarize the observations

### Data:

* '*.ast' All astronomical observations have names based on date and end with .ast
* '*.hot' Calibration requires observations of the ground with assumed temperature of 285 K
* '*.kel' Calibrated, average observations have extensions based on units.

### Directories:
* data  Selection of data for testing plotting functions.
* images Directory containing images for documenting the useage

### Support functions

The programs depend on several helper python functions:
    
* radioastronomy.py   Python to read and write spectra.  This function is shared with the data collecting software.
* interpolate.py      Python to interpolate over expected Interfering radio lines.  Needed to for more acurate calibration.
* hotcold.py	    Python to calibrate hot/cold load observations and accumulate averages
* angles.py	    Python to process angles

Examples
========

These functions must be executed in the current directry or the python programs copied to the appropriate place in your path.   

To plot all the raw data in a directory type:

R data/*

only a maximum of 25 spectra will be plotted. To plot all the hot load data

R data/*.hot


The scripts monitor the telescope azimuth and elevation and stop averaging each time the angles or
frequencies of observations change.   

The main calibration program is C

All these programs provide minimal help if executed without arguments.  Ie C

C: Calibrate Science Aficonado (NSF) horn observations
Usage: C [options]  <average_seconds> <files>

Where many parameters are optional:
-B Subtract a linear baseline fit to Spectra at Min and Max Velocities
   Min and max default velocities are:  -550.0,   180.0 km/sec
-C Flag the Center channel, using interpolation.
   This removes a strong narrow feature created by many Software Defined Radios (SDRs)
-D Additional Debug printing.
-I Optionally Flag Known Radio Frequency Interference (RFI)
   Note you need to update the c.py program to add your list of Frequencies
   RFI frequencies are Location Dependent
-H <hot load Temperature> Set the effective temperature of the hot load (Kelvins)
-N Not Calibrate.  This mode is used for tests of raw spectra
-O <output directory> Set the output directory for saved files
-R <Reference Frequency> Rest Frequency (Hz) used for Doppler calculations: 1420.406 (MHz)
-S Save average spectra in files.  The Hot and Cold Load averages are saved, too.
   Average spectra have -ave added to their names
   Calibrated spectra have a .kel (for Kelvins) extension
-X Hanning smooth the hot load observation to reduce calibration noise
-T <plot title String> Label for plot
-VA <low velocity> limit to exclude for baseline fitting
-VB <high velocity> limit to exclude for baseline fitting
Where:
   <average_seconds>: Number of seconds of observations to average.
   <average_seconds> is clock time, not observing time, so 3600. gives one plot for each hour
   <files> are Horn Observation files
   <files> must include both data pointed up (.ast) and down (.hot) observations
      All .hot files are assumed to have a system temperature of   285.0 K

 -- Glen Langston (glangsto@nsf.gov), 2018 December 11

HISTORY
18DEC11 GIL Upgrade that includes C 
18MAY04 GIL Minor corrections and place in Github
18APR30 GIL T and M plotting functions
18APR20 GIL Initial version including only the raw spectra plotting

Glen Langston, National Science Foundation (GIL)
Kevin Bandura, University of West Virginia 
=============
