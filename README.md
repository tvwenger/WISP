# WISP: Wenger Interferometry Software Package

A modular radio interferometry data calibration and imaging pipeline

## Table of Contents

* [Preamble](#preamble)

* [Caveats and Contributing](#caveats-and-contributing)

* [Requirements](#requirements)

* [Installation](#installation)

* [Calibration Quick-Start](#calibration-quick-start)

* [Imaging Quick-Start](#imaging-quick-start)

* [Configuration File](#configuration-file)

* [License and Copyright](#license-and-copyright)

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
5.0.0, 5.1.2, and 5.7.2. There may be differences between these versions
of *CASA* and previous/future versions that break the functionality
of WISP.

We have used WISP to calibrate and image continuum, spectral line, and polarization data from the Australia
Telescope Compact Array (ATCA) and the Jansky Very Large Array (JVLA),
at both ~1 GHz and ~7 GHz. If you use WISP with other telescopes or frequencies,
you may need to make modifications to the code.

If you find any bugs or would like to add any new features to WISP,
please fork the repository and submit a pull request. I will gladly
accept any corrections or contributions that extend the functionality
of WISP.

## Requirements

* Version 5 of [*CASA*](https://casa.nrao.edu/)
* LaTeX installation with the available command `pdflatex`

Note: The `matplotlib` and `astropy` packages included with *CASA*
are very out of date (at least in *CASA* version 5.7.2). You may need
to upgrade these within *CASA* via

```python
import pip
pip.main(['install','pip','--upgrade'])
pip.main(['install','astropy','--upgrade'])
pip.main(['install','matplotlib','--upgrade'])
```

## Installation

Clone this repository, then run:
```bash
/path/to/casa/bin/python setup.py install
```

## Calibration Quick-Start

Here is a short tutorial for calibrating a measurement set with
WISP.

* Create and/or edit the configuration file for your project (see
  [Configuration File](#configuration-file) section below)

* Start *CASA*

   ```bash
   /path/to/casa/bin/casa
   ```

* Import WISP

   ```python
   from wisp import wisp
   ```

* Run the WISP calibration pipeline for your measurement set (`my_data.ms`) and
  project configuration file (`my_config.ini`). You must also chose an appropriate
  reference antenna (`refant`). The arguments are `shadow_tolerance`, which is the
  tolerance in meters that an antenna may be shadowed before it is flagged, 
  `quack_interval`, which is the interval in seconds to flag at the beginning of each
  scan, `antpos` if you wish to apply antenna position corrections (VLA only),
  `gaincurve` if you wish to apply gain curve corrections, `opacity` if you wish
  to apply opacity corrections, `calpol` if you wish to apply polarization calibrations,
  `calwt` if you wish to update the data weights using the calibration solutions,
  and `solint` to set the short timescale solution interval (default is `int` for solutions
  at every integration).
  
   ```python
   wisp.calibrate('my_data.ms', 'my_config.ini', 'refant', shadow_tolerance=0.0,
                  quack_interval=0.0, antpos=True, gaincurve=True, opacity=True,
                  calpol=True, solint='10s')
   ```

* Select the required calibration tasks from the menu:

   * `0. Preliminary flags (config file, etc.) and auto-flag all fields`

         This task will flag all data matching the criteria in the
         configuration file, shadowed antennas, and the beginning and
         ends of each scan.  It will interpolate through any channels
         defined in the configuration file. Finally, it will run the
         TFCROP automatic flagging algorithm on the raw data for all
         fields.

   * `1. Auto-flag calibrator fields`

         This task will run the RFLAG automatic flagging algorithm on
         the calibrated data for all calibrator fields.

   * `2. Generate plotms figures for calibrator fields`

         This task will generate diagnostic plots for every calibrator
         field. The individual figures are saved in
         `calibrator_figures/` and compiled in
         `calibrator_figures.pdf`.

   * `3. Manually flag calibrator fields`

         This task will allow the user to interactively generate the
         PLOTMS figures and manually flag data in the calibrator
         fields. To generate an interactive PLOTMS figure from
         `calibrator_figures.pdf`, enter the Plot ID number at the
         prompt.

   * `4. Calculate and apply calibration solutions to calibrator fields`

         This task will compute all calibration tables (flux,
         bandpass, delays, gain) and apply those calibration solutions
         to the calibrator fields. The calibration solution figures
         are saved in `plotcal_figures` and compiled in
         `plotcal_figures.pdf`.

   * `5. Apply calibration solutions to science fields`

         This task takes the calibration tables from the previous task
         and applies them to the science (i.e. non-calibrator) fields.

   * `6. Auto-flag science fields`

         This task will run the RFLAG automatic flagging algorithm on
         the calibrated data for all science fields.

   * `7. Generate plotms figures for science fields`

         This task will generate diagnostic plots for every science
         field.  The individual figures are saved in
         `science_figures/` and compiled in `science_figures.pdf`.

   * `8. Manually flag science fields`

         This task will allow the user to interactively generate the
         PLOTMS figures and manually flag data in the science
         fields. To generate an interactive PLOTMS figure from
         `science_figures.pdf`, enter the Plot ID number at the
         prompt.

   * `9. Split calibrated fields`

         This task will split each calibrated field into a separated
         measurement set named `<field_name>.ms`

* Or, pass the `auto` keyword to a string of comma separated menu
  numbers to automatically run those tasks:
  
   ```python
   wisp.calibrate('my_data.ms', 'my_config.ini', 'refant', shadow_tolerance=0.0,
                  quack_interval=0.0, antpos=True, gaincurve=True, opacity=True,
                  calpol=True, solint='10s', auto='0,4,1,4,2,3')
   ```

* Here is a typical recipe for calibrating a measurement set:
   * `0.` Preliminary flagging
   * `4.` Compute calibration solutions and apply to calibrators
   * `1.` Auto-flag calibrators
   * `4.` Re-compute calibration solutions
   * `2.` Generate diagnostic plots for calibrators
   * `3.` Manually-flag calibrators
   * Repeat 4, 2, 3 until all bad data are flagged and calibration
     solutions look good.
   * `5.` Calibrate science fields
   * `6.` Auto-flag science fields
   * `7.` Generate diagnostic plots for science fields
   * `8.` Manually flag science fields
   * `9.` Split all fields

## Imaging Quick-Start

Here is a short tutorial for imaging a measurement set with WISP.

* Create and/or edit the configuration file for your project (see
  [Configuration File](#configuration-file) section below)

* Start *CASA*

   ```bash
   /path/to/casa
   ```

* Import WISP

   ```python
   from WISP import WISP
   ```

* If you are imaging a measurement set of a split field (i.e. like
  those generated by the calibration pipeline), then it is convenient to
  define a `field` variable like so. Run the WISP imaging pipeline for
  your measurement set, field, and project configuration file (`my_confit.ini`).
  The arguments are `outdir`, which is the directory in which all images
  are placed, `stokes`, which are the Stokes parameters you wish to image,
  `uvtaper` if you wish to taper the data before imaging, `outertaper`, which
  defines the size of the taper, `interactive` if you want to interactively clean,
  and `parallel` if you run to run parallel TCLEAN (note that CASA should be
  started with `mpi-casa` in this case.) 

   ```python
   field='myfieldname'
   wisp.imaging('{0}.ms'.format(field), field, 'my_config.ini', outdir='images',
                stokes='IQUV', uvtaper=True, outertaper='10arcsec', interactive=False,
                parallel=False)
   ```

* Select the required imaging tasks from the menu:

   * `0. Dirty image combined continuum spws (MFS; multi-term; multi-scale)`

         Generate a dirty image of the combined continuum spectral
         windows using multi-frequency synthesis, multi-term, and
         multi-scale TCLEAN.

   * `1. Clean combined continuum spws (MFS; multi-term; multi-scale)`

         Use auto-multithresh to automatically clean the above image.

   * `2. Dirty image each continuum spw (MFS; multi-scale)`

         Generate a multi-frequency synthesis, multi-scale TCLEAN
         dirty image of each individual continuum spectral window.

   * `3. Clean each continuum spw (MFS; multi-scale)`

         Use auto-multithresh to automatically clean the above images.

   * `4. Dirty image each continuum spw (channel; multi-scale)`

         Generate a multi-scale TCLEAN dirty data cube for each
         continuum spectral window.

   * `5. Clean each continuum spw (channel; multi-scale)`

         Use auto-multithresh to automatically clean the above cubes.

   * `6. Dirty image each line spw (MFS; multi-scale)`

         Generate a multi-frequency synthesis, multi-scale TCLEAN
         dirty image of each individual line spectral window.
      
   * `7. Clean each line spw (MFS; multi-scale)`

         Use auto-multithresh to automatically clean the above images.

   * `8. Dirty image each line spw (channel; multi-scale)`

         Generate a multi-scale TCLEAN dirty data cube for each line
         spectral window.

   * `9. Clean each line spw (channel; multi-scale)`

         Use auto-multithresh to automatically clean the above cubes.

   * `10. Generate continuum diagnostic plots`

         Generate figures for each continuum image,
         then compile those figures into `<field>.contplots.pdf`

   * `11. Generate spectral line diagnostic plots`

         Generate figures for each spectral line image,
         then compile those figures into `<field>.lineplots.pdf`.

* Or, pass the `auto` keyword to a string of comma separated menu
   numbers to automatically run those tasks:

   ```python
   field='myfieldname'
   wisp.imaging('{0}.ms'.format(field), field, 'my_config.ini', outdir='images',
                stokes='IQUV', uvtaper=True, outertaper='10arcsec', interactive=False,
                parallel=False, auto='0,1,2,3,10`)
   ```

## Configuration File

The configuration files contain the calibration and imaging parameters
that define how the data from a specific project should be processed.
A blank example configuration file is given in `config/example.ini`.
Here we explain each of the available parameters.

* [Calibrators]

   The parameters under this heading define the primary (bandpass),
   secondary (gain), flux, and polarization calibrators. If these parameters are
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

   * PolLeakage Calibrators
      
      The polarization leakage calibrator(s).

   * PolAngle Calibrators

      The polarization angle calibrator(s).

* [Calibrator Models]

   If CASA is missing the flux model (or has an incorrect model) for
   one or more of the flux calibrators, use these parameters to define
   the flux model. The calibrator models must be specified here for
   polarized calibrators. The model for flux `S` at frequency `f` is given by:
   `S(f) = S0 * (f/f0) ** (a0 + a1*log10(f/f0) + a2*log10((f/f0)**2))`
   where `S0` is the flux density at reference frequency `f0`, and `a0`, `a1`,
   and `a2` are the spectral index coefficients. The model for polarization fraction (between 0.0 and 1.0)
   is given by: `P(f) = p0 + p1 * ((f-f0)/f0) + p2 * ((f-f0)/f0)**2.0` where 
   `p0`, `p1`, and `p2` are the polarization fraction coefficients. The model for polarization angle
   (in radians) is given by: `X(f) = x0 + x1 * ((f-f0)/f0) + x2 * ((f-f0)/f0)**2.0` where 
   `x0`, `x1`, and `x2` are the polarization angle coefficients. If you wish to define
   models for multiple calibrators, do so by listing
   the values on new lines. 
   
   * Name

      The calibrator field name for which this flux model should be
      applied.

   * Refrence Frequency

      The reference frequency (with units).

   * Flux Density

      The flux density (in Jy) at the reference frequency.

   * Spectral Index Coefficients

      Comma-separated spectral index coefficients (a0,a1,a2) as
      defined in the model equation.

   * Polarization Fraction Coefficients

      Comma-separated polarization fraction coefficients (p0,p1,p2) as
      defined in the model equation.

   * Polarization Angle Coefficients

      Comma-separated polarization fraction coefficients (x0,x1,x2) as
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


* [Flags]

   These parameters define what data should be flagged throughout
   the entire measurement set. Multiple values are separated by commas,
   and the syntax follows the normal CASA FLAGDATA syntax. For example:
   ```bash
   Scan               = 0,1
   Antenna            = ea01,ea02
   Spectral Window    = 10
   Line Channels      = 0~200,900~1100
   Continuum Channels = 0,10,25,40
   ```

* [Interpolate]

   Define channels containing known bad data which should be
   overwritten with interpolated data. Values should be comma-separated,
   for example:
   
   ```bash
   Line Channels      = 256,512
   Continuum Channels = 2,5,10
   ```

* [Clean]

   Here are the various parameters associated with imaging the data.
   Most of these are explained in the [CASA manual for
   TCLEAN](https://casa.nrao.edu/docs/taskref/tclean-task.html).
   Below we describe only those not listed in the TCLEAN documentation.
   Multiple values should be comma-separated, for example:
   
   ```bash
   lineids   = H93a,H92a,H91a
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

   * lineids

     A name for each spectral line transition. Comma-separate
     multiple transitions. There needs to be one for each restfreq.

   * restfreqs

      The rest frequency (with units) for each line spectral window.
      Comma-separate multiple transitions. There needs to be one for
      each lineid.

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

## License and Copyright

GNU General Public License v3 (GNU GPLv3)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Copyright(C) 2018-2021 by
Trey V. Wenger; tvwenger@gmail.com