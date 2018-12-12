## analyze
### analyize: plotting tools for reviewing/averaging/calibrated spectra line observations from Science Aficionado Telescopes.   

The plotting programs expect Ascii input spectra with 
extensive headers describing the observations.

The plotting programs are compatible with spectra created using any of the GnuRadio NsfIntegrate??.grc designs
These python programs were created using the Gnu Radio Companion (GRC).   These programs
contained here are for analysis after the observations, not for data taking.

These programs are also compatible with the older NsfWatch output (watch.py)

### Programs:

* S     - Summarize a set of observations.  Function calls s.py
* C     - _Main_ program which calibrates, averages, flags and plots sets of spectra
* R     - Function to call r.py, which plots raw spectra

### Python:

* c.py  - Large python function to calibrate, average, flag and baseline fit observations.
* r.py  - Python function to read and plot the _raw_ spectra
* s.py  - Python function to read all selected spectra in a directory and summarize the observations

### Data:

* '*.ast' - All astronomical observations have names based on date and end with .ast
* '*.hot' - Calibration requires observations of the ground with assumed temperature of 285 K
* '*.kel' - Calibrated, average observations have extensions based on units.

### Directories:
* data      - Selection of data for testing plotting functions.  Small selections of 5 days of observations are provided in the _data_ directory to allow user testing.
* images    - Directory containing images for documenting the useage

The observations are summarized through the _S_ command.  Ie to summarize observations in the _data_ directory type:
```
S data/*

Count  Time    Az    El   G-Lon G-Lat  Frequency  BW   Gain    Filename
   1 05:01:00   0.0, 70.0 134.5, -2.6:  1421.25, 7.00  15.0 - data/18-11-01T050100.ast 
   8 12:08:31   0.0, 70.0 156.4, 43.2:  1421.25, 7.00  15.0 - data/18-11-01T120831.ast 
   1 05:00:38   0.0, 40.0 123.6, 25.5:  1421.25, 7.00  15.0 - data/18-11-02T050038.ast 
  10 12:09:52   0.0, 40.0 124.2, 27.9:  1421.25, 7.00  15.0 - data/18-11-02T120952.ast 
   1 05:00:59   0.0, 60.0 131.7,  7.1:  1421.25, 7.00  15.0 - data/18-11-03T050059.ast 
   8 12:08:40   0.0, 60.0 143.8, 39.7:  1421.25, 7.00  15.0 - data/18-11-03T120840.ast 
   1 05:00:56   0.0, 50.0 128.0, 16.4:  1421.25, 7.00  15.0 - data/18-11-04T050056.ast 
  27 17:09:42   0.0, 50.0 115.9, 36.8:  1421.25, 7.00  15.0 - data/18-11-05T170942.ast 
   1 17:40:20   0.0,-40.0 170.6,-37.5:  1421.25, 7.00  15.0 - data/18-11-05T174020.hot 
   5 17:49:20   0.0,-40.0 172.6,-35.9:  1421.25, 7.00  15.0 - data/18-11-05T174920.hot 
```
### Support functions

The programs depend on several helper python functions:
    
| code module |               Description    |
| ------------| --------- |
| radioastronomy.py | Python to read and write spectra.  This function is shared with the data collecting software.  |
| interpolate.py    | Python to interpolate over expected Interfering radio lines.  Needed to for more acurate calibration. |
| hotcold.py	    | Python to calibrate hot/cold load observations and accumulate averages. |
| angles.py	        | Python to process angle sums and differences |

These modules require several python packages, including numpy, statistics and pyephem.

![Full Calibration of 5 days of Observations, for a few minutes each day](/images/C-Cal-Baseline.png)

## Examples:

These functions must be executed in the current directry or the python programs copied to the appropriate place in your path. To plot raw data in a directory type:
```
R data/*
```
![Full Calibration of 5 days of Observations, for a few minutes each day](/images/R-spectra.png)

These data are in the _data_ subdirectory
Only a maximum of 25 spectra will be plotted.

The main calibration program is C.
To create the plot of calibrated observations (shown above) type:
```
C -B -C 4000. data/*
```
These observations were made over 5 days, with a Science Aficionados Horn and
Adalm Pluto SDR.   The Gnuradio code (see http://github.com/glangsto/gr-nsf)
was run on an Odroid XU4 octa-core computer.   The observations were setup for 7 MHz
bandwidth.  The spectra looked good, but there were some issues with full data transfer.

All these programs provide minimal help if executed without arguments.  Ie:

```
% C

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
```

Glen Langston, National Science Foundation (GIL - 2018 December 12)
