"""
imaging.py - WISP Imaging Pipeline Object

Copyright(C) 2018-2021 by
Trey V. Wenger; tvwenger@gmail.com

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

Changelog:
Trey V. Wenger - April 2021 - v3.0
    Improve code readability.
"""

import os
from astropy.coordinates import SkyCoord

import __main__ as casa


class Imaging:
    """
    The Imaging object handles all imaging steps for a single field
    in a measurement set.
    """

    def __init__(
        self,
        vis,
        field,
        logger,
        config,
        outdir=".",
        uvtaper=False,
        outertaper="",
        spws="",
        uvrange="",
        stokes="I",
        phasecenter="",
        savemodel=None,
        interactive=False,
        parallel=False,
    ):
        """
        Create a new Imaging object. Get imaging parameters from
        configuration file.

        Inputs:
          vis :: string
            The measurement set
          field :: string
            The field name to image
          logger :: logging.Logger object
            The logging object we're using
          config :: a ConfigParser object
            The ConfigParser object for this project
          outdir :: string
            The directory where to save the results. This is useful if
            you plan to make several sets of images (i.e.
            self-calibration)
          uvtaper :: boolean
            If True, get UV tapering clean parameters
          outertaper :: string
            Tapering FWHM
          spws :: string
            comma-separated list of spws to clean
            if empty, clean all spws
          uvrange :: string
            Selection on UV-range
          stokes :: string
            The Stokes parameters we're imaging. e.g. 'I' or 'IQUV'
          phasecenter :: string
            If not empty, use this phase center
          savemodel :: string
            if not none, save individual MFS images of each spectral
            window to the model column of the measurement set for
            self-calibration. This can only be done with stokes='I'.
            if savemodel == 'light': save the model after lightniter
            if savemodel == 'clean': save the model after niter
          interactive :: boolean
            if True, interactively clean
          parallel :: boolean
            if True, use parallel TCLEAN

        Returns: imaging
          imaging :: imaging.Imaging object
            a new Imaging object
        """
        self.vis = vis
        self.field = field
        self.logger = logger
        self.config = config
        self.outdir = outdir
        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)
        self.uvtaper = uvtaper
        self.outertaper = outertaper
        self.uvrange = uvrange
        self.stokes = stokes
        if savemodel is not None:
            if savemodel not in ["light", "clean"]:
                raise ValueError("Invalid savemodel: {0}".format(savemodel))
        self.savemodel = savemodel
        self.interactive = interactive
        self.parallel = parallel

        # Get continuum and line spws from configuration file
        self.cont_spws = self.config.get("Spectral Windows", "Continuum").split(",")
        self.cont_spws = [cont_spw for cont_spw in self.cont_spws if cont_spw != ""]
        self.cont_chans = ["" for _ in self.cont_spws]
        self.line_spws = self.config.get("Spectral Windows", "Line").split(",")
        self.line_spws = [line_spw for line_spw in self.line_spws if line_spw != ""]
        self.line_chans = ["" for _ in self.line_spws]
        self.logger.info("Found continuum spws: {0}".format(self.cont_spws))
        self.logger.info("Found line spws: {0}".format(self.line_spws))

        # Get clean parameters from configuration file
        self.cp = {}
        self.cp["lineids"] = config.get("Clean", "lineids").split(",")
        self.cp["restfreqs"] = config.get("Clean", "restfreqs").split(",")
        self.cp["imsize"] = [
            int(foo) for foo in config.get("Clean", "imsize").split(",")
        ]
        self.cp["frame"] = config.get("Clean", "frame")
        self.cp["pblimit"] = config.getfloat("Clean", "pblimit")
        self.cp["gridder"] = config.get("Clean", "gridder")
        self.cp["wprojplanes"] = config.getint("Clean", "wprojplanes")
        self.cp["cell"] = "{0}arcsec".format(config.getfloat("Clean", "cell"))
        self.cp["weighting"] = config.get("Clean", "weighting")
        self.cp["robust"] = config.getfloat("Clean", "robust")
        self.cp["scales"] = [
            int(foo) for foo in config.get("Clean", "scales").split(",") if foo != ""
        ]
        self.cp["gain"] = config.getfloat("Clean", "gain")
        self.cp["cyclefactor"] = config.getfloat("Clean", "cyclefactor")
        self.cp["lightniter"] = config.getint("Clean", "lightniter")
        self.cp["maxniter"] = config.getint("Clean", "maxniter")
        self.cp["nrms"] = config.getfloat("Clean", "nrms")
        self.cp["contpbchan"] = config.get("Clean", "contpbchan")
        self.cp["nterms"] = config.getint("Clean", "nterms")
        self.cp["start"] = config.get("Clean", "start")
        self.cp["width"] = config.get("Clean", "width")
        self.cp["nchan"] = config.get("Clean", "nchan")
        self.cp["chanchunks"] = config.getint("Clean", "chanchunks")
        self.cp["end"] = config.get("Clean", "end")
        self.cp["lineoutframe"] = config.get("Clean", "lineoutframe")
        self.cp["veltype"] = config.get("Clean", "veltype")
        self.cp["interpolation"] = config.get("Clean", "interpolation")
        self.cp["contoutframe"] = config.get("Clean", "contoutframe")
        self.cp["contwidth"] = config.getint("Clean", "contwidth")

        # auto-clean region parameters
        if self.uvtaper:
            heading = "Mask Taper"
        else:
            heading = "Mask NoTaper"
        self.cp["contpbmask"] = config.getfloat(heading, "contpbmask")
        self.cp["contsidelobethreshold"] = config.getfloat(
            heading, "contsidelobethreshold"
        )
        self.cp["contnoisethreshold"] = config.getfloat(heading, "contnoisethreshold")
        self.cp["contlownoisethreshold"] = config.getfloat(
            heading, "contlownoisethreshold"
        )
        self.cp["contnegativethreshold"] = config.getfloat(
            heading, "contnegativethreshold"
        )
        self.cp["contsmoothfactor"] = config.getfloat(heading, "contsmoothfactor")
        self.cp["contminbeamfrac"] = config.getfloat(heading, "contminbeamfrac")
        self.cp["contcutthreshold"] = config.getfloat(heading, "contcutthreshold")
        self.cp["contgrowiterations"] = config.getint(heading, "contgrowiterations")
        self.cp["linepbmask"] = config.getfloat(heading, "linepbmask")
        self.cp["linesidelobethreshold"] = config.getfloat(
            heading, "linesidelobethreshold"
        )
        self.cp["linenoisethreshold"] = config.getfloat(heading, "linenoisethreshold")
        self.cp["linelownoisethreshold"] = config.getfloat(
            heading, "linelownoisethreshold"
        )
        self.cp["linenegativethreshold"] = config.getfloat(
            heading, "linenegativethreshold"
        )
        self.cp["linesmoothfactor"] = config.getfloat(heading, "linesmoothfactor")
        self.cp["lineminbeamfrac"] = config.getfloat(heading, "lineminbeamfrac")
        self.cp["linecutthreshold"] = config.getfloat(heading, "linecutthreshold")
        self.cp["linegrowiterations"] = config.getint(heading, "linegrowiterations")

        # Convert phase center if necessary
        self.cp["phasecenter"] = phasecenter
        if self.cp["phasecenter"] == "" and self.cp["frame"] == "GALACTIC":
            casa.msmd.open(self.vis)
            # get field ID number
            fieldid = list(casa.msmd.fieldsforname(self.field))[0]
            phasecenter = casa.msmd.phasecenter(fieldid)
            coord = SkyCoord(
                phasecenter["m0"]["value"],
                phasecenter["m1"]["value"],
                equinox=phasecenter["refer"],
                unit=(phasecenter["m0"]["unit"], phasecenter["m1"]["unit"]),
            )
            self.cp["phasecenter"] = "GALACTIC {0}deg {1}deg".format(
                coord.galactic.l.value, coord.galactic.b.value
            )

        # Determine which spectral windows we're imaging
        if spws != "":
            my_cont_spws = []
            my_cont_chans = []
            my_line_spws = []
            my_line_chans = []
            for spwchan in spws.split(","):
                parts = spwchan.split(":")
                spw = parts[0]
                chan = ""
                if len(parts) > 1:
                    chan = parts[1]
                if spw in self.cont_spws:
                    my_cont_spws.append(spw)
                    my_cont_chans.append(chan)
                elif spw in self.line_spws:
                    my_line_spws.append(spw)
                    my_line_chans.append(chan)
                else:
                    logger.critical(
                        "Spectral window {0} not in config file".format(spw)
                    )
            self.cont_spws = my_cont_spws
            self.cont_chans = my_cont_chans
            self.line_spws = my_line_spws
            self.line_chans = my_line_chans

        # Determine which spws actually have data
        casa.msmd.open(self.vis)
        good_spws = casa.msmd.spwsforfield(self.field)
        casa.msmd.close()
        if self.cont_spws:
            self.logger.info("Checking cont spws...")
            self.cont_spws = [s for s in self.cont_spws if int(s) in good_spws]
            self.logger.info("Using cont spws: {0}".format(self.cont_spws))
        if self.line_spws:
            self.logger.info("Checking line spws...")
            self.line_spws = [s for s in self.line_spws if int(s) in good_spws]
            self.logger.info("Using line spws: {0}".format(self.line_spws))
