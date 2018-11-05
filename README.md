# WISP: Wenger Interferometry Software Package

A modular radio interferometry data calibration and imaging pipeline

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
A blank example configuration file is given in `config/example.ini`.
Here we explain each of the available parameters.

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

* [Spectral Windows]

   Here we define which spectral windows should be used for
   spectral line analyses and which should be used for continuum
   analyses. The spectral windows are comma-separated. For example:
   ```bash
   Line      = 2,5,8,11,14,17,21,22
   Continuum = 0,1,3,4,6,7,9,10,12,13,15,16,18,19,20,23
   ```

   * Line

      Which spectral windows are for spectral line analyses.

   * Continuum

      Which spectral windows are for continuum analyses.

* [Polarization]

   * Polarization

      Define the polarization direction for these data. For example:
      ```bash
      Polarization = LL,RR
      ```

* [Flags]

   These parameters define what data should be flagged throughout
   the entire measurement set. Multiple values are separated by commas,
   and the syntax follows the normal CASA FLAGDATA syntax. For example:
   ```bash
   Antenna            = ea01,ea02
   Line Channels      = 0~200,900~1100
   Continuum Channels = 0,10,25,40
   ```

   * Antenna

      Flag these antennas in all of the data.

   * Line Channels

      Flag these channels in every line spectral window.

   * Continuum Channels

      Flag these channels in every continuum spectral window.

* [Interpolate]

   Define channels containing known bad data which should be
   overwritten with interpolated data. Values should be comma-separated,
   for example:
   ```bash
   Line Channels      = 256,512
   Continuum Channels = 2,5,10
   ```

   * Line Channels

      Line spectral windows channels to interpolate.

   * Continuum Channels

      Continuum spectral window channels to interpolate.

* [Bandpass Channel Average]

   Define the number of channels to be averaged when computing
   bandpass calibration solutions. An empty value means use no
   channel averaging. For example:
   ```bash
   Line Channels      = 16
   Continuum Channels =
   ```

   * Line Channels

      The number of line spectral window channels to average.

   * Continuum Channels

      The number of continuum spectral window channels to average.

* [Clean]

   Here are the various parameters associated with imaging the data.
   Most of these are explained in the [CASA manual for
   TCLEAN](https://casa.nrao.edu/docs/taskref/tclean-task.html).
   Below we describe only those not listed in the TCLEAN documentation.
   Multiple values should be comma-separated, for example:
   ```bash
   restfreqs = 8045.605MHz,8309.385MHz,8584.823MHz
   imsize    = 600,600
   pblimit   = 0.1
   etc.
   ```

   * lightniter

      The number of CLEAN iterations used before estimating the
      image noise level.

   * maxniter

      The maximum number of CLEAN iterations to be used.

   * nrms

      The CLEAN threshold as a multiplicative factor of the
      lightly-cleaned image RMS. For example, if nrms = 2, CLEAN
      until the maximum residual is less than 2 times the
      lightly-cleaned image RMS.

   * contpbchan

      The continuum spectral window channel to use for wideband
      primary beam correction.

   * restfreqs

      The rest frequency (with units) for each line spectral window.

   * chanbuffer

      The number of channels to image beyond the range specified
      range (using start, width, nchan) to minimize the visible
      effects of interpolation. These channels are removed from the
      final image.

   * lineoutframe

      The output velocity frame for line spectral window cubes.

   * contoutframe

      The output velocity frame for continuum spectral window cubes.

* [Mask NoTaper] and [Mask Taper]

   These parameters define the automatic CLEAN-mask algorithm
   "auto-multithresh" separated for line spectral windows,
   continuum spectral windows, without uv-tapering, and with
   uv-tapering. For more information about each parameter, see
   the [auto-multithresh documentation](https://casaguides.nrao.edu/index.php/Automasking_Guide).

* [Unflag]

   These parameters define how line spectral window data should be
   automatically un-flagged when the automatic flagging algorithms
   unintentionally flag bright spectral lines. Multiple values
   are comma-separated, for example:
   ```bash
   offset = 20,20,15,15,15,10,5,0
   width  = 100
   ```
   If a spectral line is located at channel X in the last line
   spectral window, these parameters define where the spectral line
   is in every other spectral window. In each other spectral window,
   the line should be at channel X + Y where Y is the offset. The
   channels X+Y-width/2 to X+Y+width/2 are un-flagged.

   * offset

      Comma-separated channel offsets for each line spectral window
      relative to the last line spectral window. The channel offset is
      the location of a given velocity relative to the last line
      spectral window. For example, if a velocity of 0 km/s is located
      at channel 5, 10, 15, 20, and 25 in line spectral windows 1, 2,
      3, 4, and 5, the channel offsets would be -20, -15, -10, -5, and
      0.

   * width

      The total number of channels centered on the given channel to
      be un-flagged.

## Calibration Quick-Start

Here is a short tutorial for calibrating a measurement set with
WISP.

* Create and/or edit the configuration file for your project (see
  [Configuration File](#configuration-file) section below)

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
