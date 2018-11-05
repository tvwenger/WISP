# WISP: Wenger Interferometry Software Package

A modular data calibration and imaging pipeline

## Preamble

WISP is a radio interferometry calibration and imaging pipeline
written in *Python* and implemented through the *Common Astronomical
Software Applications* (CASA) Package (McMullin et al. 2007). Its
generic and modular framework is designed to handle *any* continuum or
spectral line radio interferometry data.

The main benefits of WISP over existing data reduction pipelines are
1. All of the relevant *CASA* tasks are scripted so that the user does
   not have to enter them one by one
1. WISP includes functionality to interpolate through known bad
   channels (i.e. "birdies") instead of flagging them
1. Both the calibration and imaging pipelines generate several diagnostic
   plots which are used to inspect the quality of the processing

## Caveats and Contributing

WISP has been extensively tested using *CASA* versions 4.5.2,
5.0.0, and 5.1.2. There may be differences between these versions
of *CASA* and previous/future versions that break the functionality
of WISP.

We have used WISP to calibrate and image data from the Australia
Telescope Compact Array (ATCA) and the Jansky Very Large Array (JVLA),
both a ~7 GHz. If you use WISP with other telescopes or frequencies,
you may need to make modifications to the code.

If you find any bugs or would like to add any new features to WISP,
please fork the repository and submit a pull request. I will gladly
accept any corrections or contributions that extend the functionality
of WISP.

## Requirements

* Latest version of [*CASA*](https://casa.nrao.edu/)

Note: The `matplotlib` and `astropy` packages included with *CASA*
are very out of date (at least in *CASA* version 5.1.2). You may need
to upgrade these within *CASA* via

```python
import pip
pip.main(['install','pip','--upgrade'])
pip.main(['install','astropy','--upgrade'])
pip.main(['install','matplotlib','--upgrade'])
```

## Configuration File

The configuration files contain the calibration and imaging parameters
that define how the data from a specific project should be processed.
An example calibration file is given in `config/example.ini`.  Here we
explain each of the available parameters.

* [Calibrators]

   The parameters under this heading define the primary (bandpass),
   secondary (gain), and flux calibrators. If these parameters are
   left blank, the calibrators are determined from the CASA LISTOBS
   output. If there are multiple values for a parameter, they must
   be listed on new lines. For example
   ```bash
   Primary Calibrators = 3C48
                         3C286
   ```
   * Primary Calibrators

      The primary (bandpass) calibrator(s).

   * Secondary Calibrators

      The secondary (gain) calibrator(s).

   * Flux Calibrators

      The flux calibrator(s).

* [Flux Calibrator Models]

   If CASA is missing the flux model (or has an incorrect model) for
   one or more of the flux calibrators, use these parameters to define
   the flux model. The model for flux S at frequency f is given by:
   S(f) = S0 * (f/f0) ** (a0 + a1*log10(f/f0) + a2*log10((f/f0)**2))
   where S0 is the flux density at reference frequency f0, and a0, a1,
   and a2 are the spectral index coefficients. If you wish to define
   flux models for multiple flux calibrators, do so by listing
   the values on new lines. For example
   ```bash
   Name                        = 3C48
                                 3C286
   Reference Frequency         = 5000MHz
                                 8000MHz
   Log Flux Density            = -30.8
                                 -50.2
   Spectral Index Coefficients = 26.3,-8.2,1.5
                                 40.2,-10.7,0.9
   ```
   * Name

      The calibrator field name for which this flux model should be
      applied.

   * Refrence Frequency

      The reference frequency (with units as shown in the example).

   * Log Flux Density

      The log10 of the flux density (in Jy) at the reference
      frequency.

   * Spectral Index Coefficients

      Comma-separated spectral index coefficients (a0,a1,a2) as
      defined in the model equation.

## Calibration Quick-Start

Here is a short tutorial for calibrating a measurement set with
WISP.

* Create and/or edit the configuration file for your project (see above)

* In the directory containing the measurement set, link (or copy) the
  required WISP programs
```bash
ln -s /path/to/WISP/calibration.py .
ln -s /path/to/WISP/logging.conf .
ln -s /path/to/WISP/config/my_project.ini .
```

* Start *CASA*
```bash
/path/to/casa
```

* Import the WISP calibration pipeline
```python
import calibration
```

* Run the WISP calibration pipeline for your measurement set and
  project configuration file
```python
calibration.main(vis='my_data.ms', config_file='my_project.ini')
```

* Select the required calibration tasks from the menu, or automatically
  run these steps via
```python
calibration.main(vis='my_data.ms', config_file='my_project.ini',
                 auto='0,4,1,4,5,6,7')

## Imaging Quick-Start



## Calibration Pipeline
Here are specific details about the calibration pipeline

## Imaging Pipeline
Here are specific details about the imaging pipeline
