"""
imaging.py - WISP Imaging Pipeline

Imaging a single field in a measurement set by scripting CASA tasks
and generating diagnostic plots.

Copyright(C) 2018-2019 by
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
Trey V. Wenger November 2018 - V1.0

Trey V. Wenger September 2019 - V2.0
    Added polarization imaging. Re-designed to OOP framework.
"""

import os
import glob
import re

import ConfigParser
import logging
import logging.config

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from astropy.io import fits
from astropy.wcs import WCS

import __main__ as casa

__version__ = "2.0"

# load logging configuration file
logging.config.fileConfig('logging.conf')

# catch re-naming of raw_input to input in Python3
try:
    input = raw_input
except NameError:
    pass

def natural_sort(mylist):
    """
    Natural sort an alphanumeric list

    Inputs: mylist
      mylist :: list of strings
        The list of strings to be sorted

    Returns: sorted_list
      sorted_list :: list of strings
        The sorted list
    """
    # Convert text to integers or lowercase
    convert = \
        lambda text: int(text) if text.isdigit() else text.lower()
    # define the sorting algorithm
    alphanum_key = lambda key: [convert(c) for c in
                                re.split('([0-9]+)', key)]
    # return the sorted list
    return sorted(mylist, key=alphanum_key)

class Imaging:
    """
    The Imaging object handles all imaging steps for a single field
    in a measurement set.
    """

    def __init__(self, vis, field, logger, config, uvtaper=False,
                 spws='', stokes='I', savemodel=False,
                 interactive=False):
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
          uvtaper :: boolean
            If True, get UV tapering clean parameters
          spws :: string
            comma-separated list of spws to clean
            if empty, clean all spws
          stokes :: string
            The Stokes parameters we're imaging. e.g. 'I' or 'IQUV'
          savemodel :: boolean
            if True, save individual MFS images of each spectral window
            to the model column of the measurement set for
            self-calibration. This can only be done with stokes='I'
          interactive :: boolean
            if True, interactively clean

        Returns: imaging
          imaging :: imaging.Imaging object
            a new Imaging object
        """
        self.vis = vis
        self.field = field
        self.logger = logger
        self.config = config
        self.uvtaper = uvtaper
        self.stokes = stokes
        self.savemodel = savemodel
        self.interactive = interactive
        #
        # Get continuum and line spws from configuration file
        #
        self.cont_spws = self.config.get("Spectral Windows", "Continuum")
        self.line_spws = self.config.get("Spectral Windows", "Line")
        self.logger.info("Found continuum spws: {0}".format(self.cont_spws))
        self.logger.info("Found line spws: {0}".format(self.line_spws))
        #
        # Get clean parameters from configuration file
        #
        self.cp = {}
        self.cp["lineids"] = config.get("Clean", "lineids").split(',')
        self.cp["restfreqs"] = config.get("Clean", "restfreqs").split(',')
        # clean info
        self.cp["imsize"] = [int(foo) for foo in config.get("Clean","imsize").split(',')]
        self.cp["pblimit"] = config.getfloat("Clean", "pblimit")
        self.cp["cell"] = "{0}arcsec".format(config.getfloat("Clean", "cell"))
        self.cp["weighting"] = config.get("Clean", "weighting")
        self.cp["robust"] = config.getfloat("Clean", "robust")
        self.cp["scales"] = [int(foo) for foo in config.get("Clean", "scales").split(',') if foo != '']
        self.cp["gain"] = config.getfloat("Clean", "gain")
        self.cp["cyclefactor"] = config.getfloat("Clean", "cyclefactor")
        self.cp["lightniter"] = config.getint("Clean", "lightniter")
        self.cp["maxniter"] = config.getint("Clean", "maxniter")
        self.cp["nrms"] = config.getfloat("Clean", "nrms")
        self.cp["contpbchan"] = config.get("Clean", "contpbchan")
        self.cp["nterms"] = config.getint("Clean", "nterms")
        if uvtaper:
            self.cp["outertaper"] = ["{0}arcsec".format(config.getfloat("Clean", "outertaper"))]
        else:
            self.cp["outertaper"] = []
        # re-grid info
        velstart = config.getfloat("Clean", "start")
        self.cp["velstart"] = "{0}km/s".format(velstart)
        chanwidth = config.getfloat("Clean", "width")
        self.cp["chanwidth"] = "{0}km/s".format(chanwidth)
        self.cp["nchan"] = config.getint("Clean", "nchan")
        self.cp["chanbuffer"] = config.getint("Clean", "chanbuffer")
        self.cp["cvelstart"] = "{0}km/s".format(velstart-(self.cp["chanbuffer"]*chanwidth))
        self.cp["cvelnchan"] = self.cp["nchan"]+2*self.cp["chanbuffer"]
        self.cp["lineoutframe"] = config.get("Clean", "lineoutframe")
        self.cp["veltype"] = config.get("Clean", "veltype")
        self.cp["interpolation"] = config.get("Clean", "interpolation")
        self.cp["contoutframe"] = config.get("Clean", "contoutframe")
        if uvtaper:
            heading = "Mask Taper"
        else:
            heading = "Mask NoTaper"
        self.cp["contpbmask"] = config.getfloat(heading, "contpbmask")
        self.cp["contsidelobethreshold"] = config.getfloat(heading, "contsidelobethreshold")
        self.cp["contnoisethreshold"] = config.getfloat(heading, "contnoisethreshold")
        self.cp["contlownoisethreshold"] = config.getfloat(heading, "contlownoisethreshold")
        self.cp["contnegativethreshold"] = config.getfloat(heading, "contnegativethreshold")
        self.cp["contsmoothfactor"] = config.getfloat(heading, "contsmoothfactor")
        self.cp["contminbeamfrac"] = config.getfloat(heading, "contminbeamfrac")
        self.cp["contcutthreshold"] = config.getfloat(heading, "contcutthreshold")
        self.cp["contgrowiterations"] = config.getint(heading, "contgrowiterations")
        self.cp["linepbmask"] = config.getfloat(heading, "linepbmask")
        self.cp["linesidelobethreshold"] = config.getfloat(heading, "linesidelobethreshold")
        self.cp["linenoisethreshold"] = config.getfloat(heading, "linenoisethreshold")
        self.cp["linelownoisethreshold"] = config.getfloat(heading, "linelownoisethreshold")
        self.cp["linenegativethreshold"] = config.getfloat(heading, "linenegativethreshold")
        self.cp["linesmoothfactor"] = config.getfloat(heading, "linesmoothfactor")
        self.cp["lineminbeamfrac"] = config.getfloat(heading, "lineminbeamfrac")
        self.cp["linecutthreshold"] = config.getfloat(heading, "linecutthreshold")
        self.cp["linegrowiterations"] = config.getint(heading, "linegrowiterations")
        #
        # Determine which spectral windows we're imaging
        #
        if spws != '':
            my_cont_spws = []
            my_line_spws = []
            for spw in spws.split(','):
                if spw in self.cont_spws.split(','):
                    my_cont_spws.append(spw)
                elif spw in self.line_spws.split(','):
                    my_line_spws.append(spw)
                else:
                    logger.critical('Spectral window {0} not in config file'.format(spw))
            self.cont_spws = ','.join(my_cont_spws)
            self.line_spws = ','.join(my_line_spws)
        #
        # Determine which spws actually have data
        #
        logger.info("Checking cont spws...")
        good_cont_spws = []
        for spw in self.cont_spws.split(','):
            foo = None
            foo = casa.visstat(vis=vis,spw=spw)
            if foo is not None:
                good_cont_spws.append(spw)
        self.cont_spws = ','.join(good_cont_spws)
        logger.info("Using cont spws: {0}".format(my_cont_spws))
        logger.info("Checking line spws...")
        good_line_spws = []
        for spw in self.line_spws.split(','):
            foo = None
            foo = casa.visstat(vis=vis,spw=spw)
            if foo is not None:
                good_line_spws.append(spw)
        self.line_spws = ','.join(good_line_spws)
        logger.info("Using line spws: {0}".format(my_line_spws))

    def mfs_dirty_cont(self):
        """
        Dirty image combined spws using multi-frequency synthesis

        Inputs: Nothing

        Returns: Nothing
        """
        #
        # Dirty image continuum
        #
        imagename = '{0}.cont.mfs.dirty'.format(self.field)
        if self.uvtaper:
            imagename = imagename + '.uvtaper'
        self.logger.info("Generating dirty continuum image (MFS)...")
        casa.tclean(vis=self.vis, imagename=imagename,
                    field=self.field, spw=self.cont_spws,
                    specmode='mfs', threshold='0mJy', niter=0,
                    nterms=self.cp['nterms'], deconvolver='mtmfs',
                    scales=self.cp['scales'], gain=self.cp['gain'],
                    cyclefactor=self.cp['cyclefactor'],
                    imsize=self.cp['imsize'], pblimit=-1.0,
                    cell=self.cp['cell'],
                    weighting=self.cp['weighting'],
                    robust=self.cp['robust'],
                    uvtaper=self.cp['outertaper'], stokes=self.stokes,
                    pbcor=False)
        self.logger.info("Done.")
        #
        # Primary beam correction using PB of center channel
        #
        self.logger.info("Performing primary beam correction...")
        spwlist = [int(spw) for spw in self.cont_spws.split(',')]
        weightlist = [1.0 for spw in spwlist]
        chanlist = [int(self.cp['contpbchan']) for _ in spwlist]
        casa.widebandpbcor(vis=self.vis, imagename=imagename,
                           nterms=self.cp['nterms'], pbmin=self.cp['pblimit'],
                           threshold='0.1mJy', spwlist=spwlist,
                           weightlist=weightlist, chanlist=chanlist)
        self.logger.info("Done.")
        #
        # Export to fits
        #
        self.logger.info("Exporting fits file...")
        casa.exportfits(imagename='{0}.pbcor.workdirectory/{0}.pb.tt0'.format(imagename),
                        fitsimage='{0}.cont.mfs.pb.fits'.format(self.field),
                        overwrite=True, history=False)
        casa.exportfits(imagename='{0}.image.tt0'.format(imagename),
                        fitsimage='{0}.image.fits'.format(imagename),
                        overwrite=True, history=False)
        casa.exportfits(imagename='{0}.residual.tt0'.format(imagename),
                        fitsimage='{0}.residual.fits'.format(imagename),
                        overwrite=True, history=False)
        casa.exportfits(imagename='{0}.pbcor.image.tt0'.format(imagename),
                        fitsimage='{0}.pbcor.image.fits'.format(imagename),
                        overwrite=True, history=False)
        casa.exportfits(imagename='{0}.pbcor.image.alpha'.format(imagename),
                        fitsimage='{0}.pbcor.image.alpha.fits'.format(imagename),
                        overwrite=True, history=False)
        casa.exportfits(imagename='{0}.pbcor.image.alpha.error'.format(imagename),
                        fitsimage='{0}.pbcor.image.alpha.error.fits'.format(imagename),
                        overwrite=True, history=False)
        self.logger.info("Done.")

    def mfs_clean_cont(self):
        """
        Clean combined spws using multi-frequency synthesis

        Inputs: Nothing

        Returns: Nothing
        """
        #
        # If not cleaning interactively, lightly clean to get RMS
        # threshold
        #
        if not self.interactive:
            imagename = '{0}.cont.mfs.lightclean'.format(self.field)
            if self.uvtaper:
                imagename = imagename + '.uvtaper'
            self.logger.info("Lightly cleaning continuum image (MFS)...")
            casa.tclean(vis=self.vis, imagename=imagename,
                        field=self.field, spw=self.cont_spws,
                        specmode='mfs', threshold='0mJy',
                        niter=self.cp["lightniter"],
                        usemask='auto-multithresh',
                        pbmask=self.cp['contpbmask'],
                        sidelobethreshold=self.cp['contsidelobethreshold'],
                        noisethreshold=self.cp['contnoisethreshold'],
                        lownoisethreshold=self.cp['contlownoisethreshold'],
                        negativethreshold=self.cp['contnegativethreshold'],
                        smoothfactor=self.cp['contsmoothfactor'],
                        minbeamfrac=self.cp['contminbeamfrac'],
                        cutthreshold=self.cp['contcutthreshold'],
                        growiterations=self.cp['contgrowiterations'],
                        nterms=self.cp['nterms'], deconvolver='mtmfs',
                        scales=self.cp['scales'],
                        gain=self.cp['gain'],
                        cyclefactor=self.cp['cyclefactor'],
                        imsize=self.cp['imsize'], pblimit=-1.0,
                        cell=self.cp['cell'],
                        weighting=self.cp['weighting'],
                        robust=self.cp['robust'],
                        uvtaper=self.cp['outertaper'],
                        stokes=self.stokes, pbcor=False)
            self.logger.info("Done.")
            #
            # Get RMS of residuals outside of clean mask
            #
            dat = casa.imstat(imagename='{0}.residual.tt0'.format(imagename),
                              axes=[0, 1], mask="'{0}.mask' == 0".format(imagename))
            self.logger.info("Max un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.max(dat['rms'])))
            self.logger.info("Mean un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.mean(dat['rms'])))
            self.logger.info("Median un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.median(dat['rms'])))
            self.logger.info("Max un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.max(dat['medabsdevmed'])))
            self.logger.info("Mean un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.mean(dat['medabsdevmed'])))
            self.logger.info("Median un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.median(dat['medabsdevmed'])))
            self.logger.info("Using median MADS*1.4826 times {0} (user-defined) as threshold".format(self.cp['nrms']))
            threshold = '{0:.2f}mJy'.format(self.cp["nrms"]*1000.*1.4826*np.median(dat['medabsdevmed']))
        else:
            #
            # No threshold for interactive clean
            #
            threshold = '0.0mJy'
        #
        # Clean to threshold
        #
        imagename = '{0}.cont.mfs.clean'.format(self.field)
        if self.uvtaper:
            imagename = imagename + '.uvtaper'
        self.logger.info("Cleaning continuum image (MFS) to threshold: {0}...".format(threshold))
        casa.tclean(vis=self.vis, imagename=imagename,
                    field=self.field, spw=self.cont_spws,
                    specmode='mfs', threshold=threshold,
                    niter=self.cp['maxniter'],
                    usemask='auto-multithresh',
                    pbmask=self.cp['contpbmask'],
                    sidelobethreshold=self.cp['contsidelobethreshold'],
                    noisethreshold=self.cp['contnoisethreshold'],
                    lownoisethreshold=self.cp['contlownoisethreshold'],
                    negativethreshold=self.cp['contnegativethreshold'],
                    smoothfactor=self.cp['contsmoothfactor'],
                    minbeamfrac=self.cp['contminbeamfrac'],
                    cutthreshold=self.cp['contcutthreshold'],
                    growiterations=self.cp['contgrowiterations'],
                    nterms=self.cp['nterms'], deconvolver='mtmfs',
                    scales=self.cp['scales'], gain=self.cp['gain'],
                    cyclefactor=self.cp['cyclefactor'],
                    imsize=self.cp['imsize'], pblimit=-1.0,
                    cell=self.cp['cell'],
                    weighting=self.cp['weighting'],
                    robust=self.cp['robust'],
                    uvtaper=self.cp['outertaper'], pbcor=False,
                    stokes=self.stokes, interactive=self.interactive)
        self.logger.info("Done.")
        #
        # Primary beam correction using PB of center channel
        #
        self.logger.info("Performing primary beam correction...")
        spwlist = [int(spw) for spw in self.cont_spws.split(',')]
        weightlist = [1.0 for spw in spwlist]
        chanlist = [int(self.cp['contpbchan']) for _ in spwlist]
        casa.widebandpbcor(vis=self.vis, imagename=imagename,
                           nterms=self.cp['nterms'],
                           pbmin=self.cp['pblimit'],
                           threshold='0.1mJy', spwlist=spwlist,
                           weightlist=weightlist, chanlist=chanlist)
        self.logger.info("Done.")
        #
        # Export to fits
        #
        self.logger.info("Exporting fits file...")
        casa.exportfits(imagename='{0}.pbcor.workdirectory/{0}.pb.tt0'.format(imagename),
                        fitsimage='{0}.pb.fits'.format(imagename),
                        overwrite=True, history=False)
        casa.exportfits(imagename='{0}.image.tt0'.format(imagename),
                        fitsimage='{0}.image.fits'.format(imagename),
                        overwrite=True, history=False)
        casa.exportfits(imagename='{0}.mask'.format(imagename),
                        fitsimage='{0}.mask.fits'.format(imagename),
                        overwrite=True, history=False)
        casa.exportfits(imagename='{0}.residual.tt0'.format(imagename),
                        fitsimage='{0}.residual.fits'.format(imagename),
                        overwrite=True, history=False)
        casa.exportfits(imagename='{0}.pbcor.image.tt0'.format(imagename),
                        fitsimage='{0}.pbcor.image.fits'.format(imagename),
                        overwrite=True, history=False)
        casa.exportfits(imagename='{0}.pbcor.image.alpha'.format(imagename),
                        fitsimage='{0}.pbcor.image.alpha.fits'.format(imagename),
                        overwrite=True, history=False)
        casa.exportfits(imagename='{0}.pbcor.image.alpha.error'.format(imagename),
                        fitsimage='{0}.pbcor.image.alpha.error.fits'.format(imagename),
                        overwrite=True, history=False)
        self.logger.info("Done.")

    def mfs_dirty_spws(self, spws):
        """
        Dirty image each supplied spw (MFS)

        Inputs:
          spws :: string
            comma-separated list of spws to image

        Returns: Nothing
        """
        for spw in spws.split(','):
            #
            # Dirty image
            #
            imagename = '{0}.spw{1}.mfs.dirty'.format(self.field, spw)
            if self.uvtaper:
                imagename = imagename + '.uvtaper'
            self.logger.info("Generating dirty image of spw {0} (MFS)...".format(spw))
            casa.tclean(vis=self.vis, imagename=imagename,
                        field=self.field, spw=spw, specmode='mfs',
                        threshold='0mJy', niter=0,
                        deconvolver='multiscale',
                        scales=self.cp['scales'],
                        gain=self.cp['gain'],
                        cyclefactor=self.cp['cyclefactor'],
                        imsize=self.cp['imsize'], pblimit=-1.0,
                        cell=self.cp['cell'],
                        weighting=self.cp['weighting'],
                        robust=self.cp['robust'],
                        uvtaper=self.cp['outertaper'],
                        stokes=self.stokes, pbcor=False)
            self.logger.info("Done.")
            #
            # Generate primary beam image
            #
            pbimagename = '{0}.spw{1}.mfs.pb'.format(self.field, spw)
            if self.uvtaper:
                pbimagename = pbimagename + '.uvtaper'
            self.logger.info("Generating primary beam image of spw {0} (MFS)...".format(spw))
            casa.tclean(vis=self.vis, imagename=pbimagename,
                        field=self.field, spw=spw, specmode='mfs',
                        threshold='0mJy', niter=0,
                        deconvolver='multiscale',
                        scales=self.cp['scales'],
                        gain=self.cp['gain'],
                        cyclefactor=self.cp['cyclefactor'],
                        imsize=self.cp['imsize'],
                        pblimit=self.cp['pblimit'],
                        cell=self.cp['cell'],
                        weighting=self.cp['weighting'],
                        robust=self.cp['robust'],
                        uvtaper=self.cp['outertaper'], pbcor=False)
            self.logger.info("Done.")
            #
            # Primary beam correction
            #
            self.logger.info("Performing primary beam correction...")
            casa.impbcor(imagename='{0}.image'.format(imagename),
                         pbimage='{0}.pb'.format(pbimagename),
                         outfile='{0}.pbcor.image'.format(imagename),
                         overwrite=True)
            self.logger.info("Done.")
            #
            # Export to fits
            #
            self.logger.info("Exporting fits file...")
            casa.exportfits(imagename='{0}.pb'.format(pbimagename),
                            fitsimage='{0}.fits'.format(pbimagename),
                            overwrite=True, history=False)
            casa.exportfits(imagename='{0}.image'.format(imagename),
                            fitsimage='{0}.image.fits'.format(imagename),
                            overwrite=True, history=False)
            casa.exportfits(imagename='{0}.residual'.format(imagename),
                            fitsimage='{0}.residual.fits'.format(imagename),
                            overwrite=True, history=False)
            casa.exportfits(imagename='{0}.pbcor.image'.format(imagename),
                            fitsimage='{0}.pbcor.image.fits'.format(imagename),
                            overwrite=True, history=False)
            self.logger.info("Done.")

    def mfs_clean_spws(self, spws):
        """
        Clean each supplied spw (MFS)

        Inputs:
          spws :: string
            comma-separated string of spws to image

        Returns: Nothing
        """
        for spw in spws.split(','):
            #
            # If not interactive, Lightly clean to get threshold
            #
            if not self.interactive:
                imagename = '{0}.spw{1}.mfs.lightclean'.format(self.field, spw)
                if self.uvtaper:
                    imagename = imagename + '.uvtaper'
                self.logger.info("Lightly cleaning spw {0} (MFS)...".format(spw))
                casa.tclean(vis=self.vis, imagename=imagename,
                            field=self.field, spw=spw, specmode='mfs',
                            threshold='0mJy',
                            niter=self.cp["lightniter"],
                            usemask='auto-multithresh',
                            pbmask=self.cp['linepbmask'],
                            sidelobethreshold=self.cp['linesidelobethreshold'],
                            noisethreshold=self.cp['linenoisethreshold'],
                            lownoisethreshold=self.cp['linelownoisethreshold'],
                            negativethreshold=self.cp['linenegativethreshold'],
                            smoothfactor=self.cp['linesmoothfactor'],
                            minbeamfrac=self.cp['lineminbeamfrac'],
                            cutthreshold=self.cp['linecutthreshold'],
                            growiterations=self.cp['linegrowiterations'],
                            deconvolver='multiscale',
                            scales=self.cp['scales'],
                            gain=self.cp['gain'],
                            cyclefactor=self.cp['cyclefactor'],
                            imsize=self.cp['imsize'], pblimit=-1.0,
                            cell=self.cp['cell'],
                            weighting=self.cp['weighting'],
                            robust=self.cp['robust'],
                            uvtaper=self.cp['outertaper'],
                            stokes=self.stokes, pbcor=False)
                self.logger.info("Done.")
                #
                # Get RMS of residuals
                #
                dat = casa.imstat(imagename='{0}.residual'.format(imagename),
                                  axes=[0, 1], mask="'{0}.mask' == 0".format(imagename))
                self.logger.info("Max un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.max(dat['rms'])))
                self.logger.info("Mean un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.mean(dat['rms'])))
                self.logger.info("Median un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.median(dat['rms'])))
                self.logger.info("Max un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.max(dat['medabsdevmed'])))
                self.logger.info("Mean un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.mean(dat['medabsdevmed'])))
                self.logger.info("Median un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.median(dat['medabsdevmed'])))
                self.logger.info("Using median MADS*1.4826 times {0} (user-defined) as threshold".format(self.cp['nrms']))
                threshold = '{0:.2f}mJy'.format(self.cp["nrms"]*1000.*1.4826*np.median(dat['medabsdevmed']))
            else:
                threshold = '0.0mJy'
            #
            # Save model if necessary
            #
            savemodel = 'none'
            if self.savemodel:
                savemodel = 'modelcolumn'
            #
            # Clean to threshold
            #
            imagename = '{0}.spw{1}.mfs.clean'.format(self.field, spw)
            if self.uvtaper:
                imagename = imagename + '.uvtaper'
            self.logger.info("Cleaning spw {0} (MFS) to threshold: {1}...".format(spw, threshold))
            casa.tclean(vis=vis, imagename=imagename, field=self.field,
                        spw=spw, specmode='mfs', threshold=threshold,
                        niter=self.cp["maxniter"],
                        usemask='auto-multithresh',
                        pbmask=self.cp['linepbmask'],
                        sidelobethreshold=self.cp['linesidelobethreshold'],
                        noisethreshold=self.cp['linenoisethreshold'],
                        lownoisethreshold=self.cp['linelownoisethreshold'],
                        negativethreshold=self.cp['linenegativethreshold'],
                        smoothfactor=self.cp['linesmoothfactor'],
                        minbeamfrac=self.cp['lineminbeamfrac'],
                        cutthreshold=self.cp['linecutthreshold'],
                        growiterations=self.cp['linegrowiterations'],
                        deconvolver='multiscale',
                        scales=self.cp['scales'],
                        gain=self.cp['gain'],
                        cyclefactor=self.cp['cyclefactor'],
                        imsize=self.cp['imsize'], pblimit=-1.0,
                        cell=self.cp['cell'],
                        weighting=self.cp['weighting'],
                        robust=self.cp['robust'],
                        uvtaper=self.cp['outertaper'], pbcor=False,
                        stokes=self.stokes, savemodel=savemodel,
                        interactive=self.interactive)
            self.logger.info("Done.")
            #
            # Primary beam correction
            #
            pbimagename = '{0}.spw{1}.mfs.pb'.format(self.field, spw)
            if self.uvtaper:
                pbimagename = pbimagename + '.uvtaper'
            self.logger.info("Performing primary beam correction...")
            casa.impbcor(imagename='{0}.image'.format(imagename),
                         pbimage='{0}.pb'.format(pbimagename),
                         outfile='{0}.pbcor.image'.format(imagename),
                         overwrite=True)
            self.logger.info("Done.")
            #
            # Export to fits
            #
            self.logger.info("Exporting fits file...")
            casa.exportfits(imagename='{0}.pb'.format(pbimagename),
                            fitsimage='{0}.fits'.format(pbimagename),
                            overwrite=True, history=False)
            casa.exportfits(imagename='{0}.image'.format(imagename),
                            fitsimage='{0}.image.fits'.format(imagename),
                            overwrite=True, history=False)
            casa.exportfits(imagename='{0}.mask'.format(imagename),
                            fitsimage='{0}.mask.fits'.format(imagename),
                            overwrite=True, history=False)
            casa.exportfits(imagename='{0}.residual'.format(imagename),
                            fitsimage='{0}.residual.fits'.format(imagename),
                            overwrite=True, history=False)
            casa.exportfits(imagename='{0}.pbcor.image'.format(imagename),
                            fitsimage='{0}.pbcor.image.fits'.format(imagename),
                            overwrite=True, history=False)
            self.logger.info("Done.")

    def channel_dirty_spws(self, spws, spwtype):
        """
        Dirty image each supplied spectral window (channel cube)

        Inputs:
          spws :: string
            comma-separated string of spws to image
          spwtype :: string
            'line' or 'cont', to determine which clean params to use

        Returns: Nothing
        """
        #
        # Set channel parameters
        #
        if spwtype == 'cont':
            restfreqs = [None for spw in spws.split(',')]
            start = None
            width = None
            nchan = None
            outframe = self.cp['contoutframe']
            veltype = None
            interpolation = None
        elif spwtype == 'line':
            spw_inds = [self.config.get("Spectral Windows","Line").split(',').index(spw)
                        for spw in spws.split(',')]
            restfreqs = [self.config.get("Clean","restfreqs").split(',')[spw_ind]
                         for spw_ind in spw_inds]
            start = self.cp['velstart']
            width = self.cp['chanwidth']
            nchan = self.cp['nchan']
            outframe = self.cp['lineoutframe']
            veltype = self.cp['veltype']
            interpolation = self.cp['interpolation']
        else:
            self.logger.critical("Error: spwtype {0} not supported".format(spwtype))
            raise ValueError("Invalid spwtype")
        #
        # Loop over spws
        #
        for spw, restfreq in zip(spws.split(','), restfreqs):
            #
            # dirty image spw
            #
            imagename = '{0}.spw{1}.channel.dirty'.format(field, spw)
            if self.uvtaper:
                imagename = imagename + '.uvtaper'
            self.logger.info("Dirty imaging spw {0} (restfreq: {1})...".format(spw,restfreq))
            casa.tclean(vis=self.vis, imagename=imagename,
                        field=self.field, spw=spw, specmode='cube',
                        threshold='0mJy', niter=0,
                        deconvolver='multiscale',
                        scales=self.cp['scales'],
                        gain=self.cp['gain'],
                        cyclefactor=self.cp['cyclefactor'],
                        imsize=self.cp['imsize'], pblimit=-1.0,
                        cell=self.cp['cell'],
                        weighting=self.cp['weighting'],
                        robust=self.cp['robust'], restfreq=restfreq,
                        start=start, width=width, nchan=nchan,
                        outframe=outframe, veltype=veltype,
                        interpolation=interpolation,
                        uvtaper=self.cp['outertaper'],
                        stokes=self.stokes, pbcor=False)
            self.logger.info("Done.")
            #
            # Generate primary beam image
            #
            pbimagename = '{0}.spw{1}.channel.pb'.format(self.field, spw)
            if self.uvtaper:
                pbimagename = pbimagename + '.uvtaper'
            self.logger.info("Generating primary beam image of spw {0} (MFS)...".format(spw))
            casa.tclean(vis=self.vis, imagename=pbimagename,
                        field=self.field, spw=spw, specmode='cube',
                        threshold='0mJy', niter=0,
                        deconvolver='multiscale',
                        scales=self.cp['scales'],
                        gain=self.cp['gain'],
                        cyclefactor=self.cp['cyclefactor'],
                        imsize=self.cp['imsize'],
                        pblimit=self.cp['pblimit'],
                        cell=self.cp['cell'],
                        weighting=self.cp['weighting'],
                        robust=self.cp['robust'], restfreq=restfreq,
                        start=start, width=width, nchan=nchan,
                        outframe=outframe, veltype=veltype,
                        interpolation=interpolation,
                        uvtaper=self.cp['outertaper'],
                        pbcor=False)
            self.logger.info("Done.")
            #
            # Primary beam correction
            #
            self.logger.info("Performing primary beam correction...")
            casa.impbcor(imagename='{0}.image'.format(imagename),
                         pbimage='{0}.pb'.format(pbimagename),
                         outfile='{0}.pbcor.image'.format(imagename),
                         overwrite=True)
            self.logger.info("Done.")
            #
            # Export to fits
            #
            self.logger.info("Exporting fits file...")
            casa.exportfits(imagename='{0}.image'.format(imagename),
                            fitsimage='{0}.image.fits'.format(imagename),
                            velocity=True, overwrite=True,
                            history=False)
            casa.exportfits(imagename='{0}.residual'.format(imagename),
                            fitsimage='{0}.residual.fits'.format(imagename),
                            velocity=True, overwrite=True,
                            history=False)
            casa.exportfits(imagename='{0}.pbcor.image'.format(imagename),
                            fitsimage='{0}.pbcor.image.fits'.format(imagename),
                            velocity=True, overwrite=True,
                            history=False)
            self.logger.info("Done.")

    def channel_clean_spws(self, spws, spwtype):
        """
        Clean all supplied spws by channel using clean mask from MFS
        images.

        Inputs:
          spws :: string
            comma-separated string of spws to image
          spwtype :: string
            'line' or 'cont', to determine which clean params to use

        Returns: Nothing
        """
        #
        # Set channel parameters
        #
        if spwtype == 'cont':
            restfreqs = [None for spw in spws.split(',')]
            start = None
            width = None
            nchan = None
            outframe = self.cp['contoutframe']
            veltype = None
            interpolation = None
        elif spwtype == 'line':
            spw_inds = [self.config.get("Spectral Windows","Line").split(',').index(spw)
                        for spw in spws.split(',')]
            restfreqs = [self.config.get("Clean","restfreqs").split(',')[spw_ind]
                         for spw_ind in spw_inds]
            start = self.cp['velstart']
            width = self.cp['chanwidth']
            nchan = self.cp['nchan']
            outframe = self.cp['lineoutframe']
            veltype = self.cp['veltype']
            interpolation = self.cp['interpolation']
        else:
            self.logger.critical("Error: spwtype {0} not supported".format(spwtype))
            raise ValueError("Invalid spwtype")
        #
        # Loop over spws
        #
        for spw, restfreq in zip(spws.split(','), restfreqs):
            #
            # If not interactive, Lightly clean spw
            #
            if not self.interactive:
                imagename = '{0}.spw{1}.channel.lightclean'.format(self.field, spw)
                if self.uvtaper:
                    imagename = imagename + '.uvtaper'
                    mask = '{0}.spw{1}.mfs.clean.uvtaper.mask'.format(self.field, spw)
                else:
                    mask = '{0}.spw{1}.mfs.clean.mask'.format(self.field, spw)
                if not os.path.isdir(mask):
                    self.logger.critical("Error: {0} does not exist".format(mask))
                    raise ValueError("{0} does not exist".format(mask))
                self.logger.info("Lightly cleaning spw {0} (restfreq: {1})...".format(spw, restfreq))
                self.logger.info("Using mask: {0}".format(mask))
                casa.tclean(vis=self.vis, imagename=imagename,
                            field=self.field, spw=spw,
                            specmode='cube', threshold='0mJy',
                            niter=self.cp['lightniter']*self.cp['nchan'],
                            mask=mask, deconvolver='multiscale',
                            scales=self.cp['scales'],
                            gain=self.cp['gain'],
                            cyclefactor=self.cp['cyclefactor'],
                            imsize=self.cp['imsize'], pblimit=-1.0,
                            cell=self.cp['cell'],
                            weighting=self.cp['weighting'],
                            robust=self.cp['robust'],
                            restfreq=restfreq, start=start,
                            width=width, nchan=nchan,
                            outframe=outframe, veltype=veltype,
                            interpolation=interpolation,
                            uvtaper=self.cp['outertaper'],
                            stokes=self.stokes, pbcor=False)
                #
                # Get RMS of residuals
                #
                dat = casa.imstat(imagename='{0}.residual'.format(imagename),
                                  axes=[0, 1], mask="'{0}.mask' == 0".format(imagename))
                self.logger.info("Max un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.max(dat['rms'])))
                self.logger.info("Mean un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.mean(dat['rms'])))
                self.logger.info("Median un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.median(dat['rms'])))
                self.logger.info("Max un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.max(dat['medabsdevmed'])))
                self.logger.info("Mean un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.mean(dat['medabsdevmed'])))
                self.logger.info("Median un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.median(dat['medabsdevmed'])))
                self.logger.info("Using median MADS*1.4826 times {0} (user-defined) as threshold".format(self.cp['nrms']))
                threshold = '{0:.2f}mJy'.format(self.cp["nrms"]*1000.*1.4826*np.median(dat['medabsdevmed']))
            else:
                threshold = '0.0mJy'
            #
            # Deep clean to threshold
            #
            imagename = '{0}.spw{1}.channel.clean'.format(self.field, spw)
            if self.uvtaper:
                imagename = imagename + '.uvtaper'
                mask = '{0}.spw{1}.mfs.clean.uvtaper.mask'.format(self.field, spw)
            else:
                mask = '{0}.spw{1}.mfs.clean.mask'.format(self.field, spw)
            if not os.path.isdir(mask):
                self.logger.critical("Error: {0} does not exist".format(mask))
                raise ValueError("{0} does not exist".format(mask))
            self.logger.info("Cleaning spw {0} (restfreq: {1}) to threshold: {2}...".format(spw, restfreq, threshold))
            self.logger.info("Using mask: {0}".format(mask))
            casa.tclean(vis=self.vis, imagename=imagename,
                        field=self.field, spw=spw, specmode='cube',
                        threshold=threshold,
                        niter=self.cp['maxniter']*self.cp['nchan'],
                        mask=mask, deconvolver='multiscale',
                        scales=self.cp['scales'],
                        gain=self.cp['gain'],
                        cyclefactor=self.cp['cyclefactor'],
                        imsize=self.cp['imsize'], pblimit=-1.0,
                        cell=self.cp['cell'],
                        weighting=self.cp['weighting'],
                        robust=self.cp['robust'], restfreq=restfreq,
                        start=start, width=width, nchan=nchan,
                        outframe=outframe, veltype=veltype,
                        interpolation=interpolation,
                        uvtaper=self.cp['outertaper'], pbcor=False,
                        stokes=self.stokes,
                        interactive=self.interactive)
            self.logger.info("Done.")
            #
            # Primary beam correction
            #
            pbimagename = '{0}.spw{1}.channel.pb'.format(self.field, spw)
            if self.uvtaper:
                pbimagename = pbimagename + '.uvtaper'
            self.logger.info("Performing primary beam correction...")
            casa.impbcor(imagename='{0}.image'.format(imagename),
                         pbimage='{0}.pb'.format(pbimagename),
                         outfile='{0}.pbcor.image'.format(imagename),
                         overwrite=True)
            self.logger.info("Done.")
            #
            # Export to fits
            #
            self.logger.info("Exporting fits file...")
            casa.exportfits(imagename='{0}.image'.format(imagename),
                            fitsimage='{0}.image.fits'.format(imagename),
                            velocity=True, overwrite=True,
                            history=False)
            casa.exportfits(imagename='{0}.mask'.format(imagename),
                            fitsimage='{0}.mask.fits'.format(imagename),
                            overwrite=True, history=False)
            casa.exportfits(imagename='{0}.residual'.format(imagename),
                            fitsimage='{0}.residual.fits'.format(imagename),
                            velocity=True, overwrite=True,
                            history=False)
            casa.exportfits(imagename='{0}.pbcor.image'.format(imagename),
                            fitsimage='{0}.pbcor.image.fits'.format(imagename),
                            velocity=True, overwrite=True,
                            history=False)
            self.logger.info("Done.")

    def contplot(self):
        """
        Generate PDF of MFS diagnostic plots

        Inputs: Nothing

        Returns: Nothing
        """
        #
        # Get center pixels
        #
        center_x = int(self.cp['imsize'][0]/2)
        center_y = int(self.cp['imsize'][1]/2)
        #
        # Loop over all plot filenames
        #
        goodplots = []
        spws = ['cont'] + natural_sort(self.cont_spws.split(',')+
                                       self.line_spws.split(','))
        for spw in spws:
            if spw != 'cont':
                spw = 'spw{0}'.format(spw)
            # check that this spectral window exists
            fname = '{0}.{1}.mfs.clean.image.fits'.format(self.field, spw)
            if self.uvtaper:
                fname = '{0}.{1}.mfs.clean.uvtaper.image.fits'.format(self.field, spw)
            if not os.path.exists(fname):
                continue
            if self.uvtaper:
                fitsfiles = ['{0}.{1}.mfs.dirty.uvtaper.image.fits'.format(self.field, spw),
                             '{0}.{1}.mfs.clean.uvtaper.image.fits'.format(self.field, spw),
                             '{0}.{1}.mfs.clean.uvtaper.residual.fits'.format(self.field, spw),
                             '{0}.{1}.mfs.clean.uvtaper.pbcor.image.fits'.format(self.field, spw)]
                maskfiles = ['{0}.{1}.mfs.clean.uvtaper.mask.fits'.format(self.field, spw),
                             '{0}.{1}.mfs.clean.uvtaper.mask.fits'.format(self.field, spw),
                             '{0}.{1}.mfs.clean.uvtaper.mask.fits'.format(self.field, spw),
                             '{0}.{1}.mfs.clean.uvtaper.mask.fits'.format(self.field, spw)]
                titles = ['{0} - spw: {1} - UVTap/Dirty'.format(self.field, spw),
                          '{0} - spw: {1} - UVTap/Clean'.format(self.field, spw),
                          '{0} - spw: {1} - UVTap/Residual'.format(self.field, spw),
                          '{0} - spw: {1} - UVTap/PBCorr'.format(self.field, spw)]
                labels = ['Flux Density (Jy/beam)',
                          'Flux Density (Jy/beam)',
                          'Flux Density (Jy/beam)',
                          'Flux Density (Jy/beam)']
                vlims = [(None, None),
                         (None, None),
                         (None, None),
                         (None, None)]
            else:
                fitsfiles = ['{0}.{1}.mfs.dirty.image.fits'.format(self.field, spw),
                             '{0}.{1}.mfs.clean.image.fits'.format(self.field, spw),
                             '{0}.{1}.mfs.clean.residual.fits'.format(self.field, spw),
                             '{0}.{1}.mfs.clean.pbcor.image.fits'.format(self.field, spw)]
                maskfiles = ['{0}.{1}.mfs.clean.mask.fits'.format(self.field, spw),
                             '{0}.{1}.mfs.clean.mask.fits'.format(self.field, spw),
                             '{0}.{1}.mfs.clean.mask.fits'.format(self.field, spw),
                             '{0}.{1}.mfs.clean.mask.fits'.format(self.field, spw)]
                titles = ['{0} - spw: {1} - Dirty'.format(self.field, spw),
                          '{0} - spw: {1} - Clean'.format(self.field, spw),
                          '{0} - spw: {1} - Residual'.format(self.field, spw),
                          '{0} - spw: {1} - PBCorr'.format(self.field, spw)]
                labels = ['Flux Density (Jy/beam)',
                          'Flux Density (Jy/beam)',
                          'Flux Density (Jy/beam)',
                          'Flux Density (Jy/beam)']
                vlims = [(None, None),
                         (None, None),
                         (None, None),
                         (None, None)]
            for fitsfile, maskfile, title, label, vlim in \
                zip(fitsfiles, maskfiles, titles, labels, vlims):
                #
                # Open fits file, generate WCS
                #
                hdu = fits.open(fitsfile)[0]
                wcs = WCS(hdu.header)
                #
                # Generate figure
                #
                plt.ioff()
                fig = plt.figure()
                ax = plt.subplot(projection=wcs.sub(['celestial']))
                ax.set_title(title)
                cax = ax.imshow(hdu.data[0, 0],
                                origin='lower', interpolation='none',
                                cmap='binary', vmin=vlim[0],
                                vmax=vlim[1])
                # ax.grid(True,color='black',ls='solid')
                ax.coords[0].set_major_formatter('hh:mm:ss')
                ax.set_xlabel('RA (J2000)')
                ax.set_ylabel('Declination (J2000)')
                #
                # Plot clean mask as contour
                #
                if maskfile != '':
                    maskhdu = fits.open(maskfile)[0]
                    contour = ax.contour(maskhdu.data[0,0],levels=[0.5],
                                         origin='lower',colors='r',linewidths=2)
                #
                # Plot beam, if it is defined
                #
                if 'BMAJ' in hdu.header.keys():
                    cell = float(self.cp['cell'].replace('arcsec', ''))
                    beam_maj = hdu.header['BMAJ']*3600./cell
                    beam_min = hdu.header['BMIN']*3600./cell
                    beam_pa = hdu.header['BPA']
                    ellipse = Ellipse((center_x-int(3.*center_x/4),
                                       center_y-int(3.*center_y/4)),
                                      beam_min, beam_maj,
                                      angle=beam_pa, fill=True,
                                      zorder=10, hatch='///',
                                      edgecolor='black',
                                      facecolor='white')
                    ax.add_patch(ellipse)
                #
                # Plot colorbar
                #
                cbar = fig.colorbar(cax, fraction=0.046, pad=0.04)
                cbar.set_label(label)
                #
                # Re-scale to fit, then save
                #
                fname = fitsfile.replace('.fits', '.pdf')
                fig.savefig(fname, bbox_inches='tight')
                plt.close(fig)
                plt.ion()
                goodplots.append(fname)
        #
        # Generate PDF of plots
        #
        # need to fix filenames so LaTeX doesn't complain
        outplots = ['{'+fn.replace('.pdf', '')+'}.pdf' for fn in goodplots]
        self.logger.info("Generating PDF...")
        fname = '{0}.contplots.tex'.format(self.field)
        if self.uvtaper:
            fname = '{0}.uvtaper.contplots.tex'.format(self.field)
        with open(fname, 'w') as f:
            f.write(r"\documentclass{article}"+"\n")
            f.write(r"\usepackage{graphicx}"+"\n")
            f.write(r"\usepackage[margin=0.1cm]{geometry}"+"\n")
            f.write(r"\begin{document}"+"\n")
            for i in range(0, len(outplots), 4):
                f.write(r"\begin{figure}"+"\n")
                f.write(r"\centering"+"\n")
                f.write(r"\includegraphics[width=0.45\textwidth]{"+outplots[i]+r"}" + "\n")
                f.write(r"\includegraphics[width=0.45\textwidth]{"+outplots[i+1]+r"} \\" + "\n")
                f.write(r"\includegraphics[width=0.45\textwidth]{"+outplots[i+2]+r"}" + "\n")
                f.write(r"\includegraphics[width=0.45\textwidth]{"+outplots[i+3]+r"}" + "\n")
                f.write(r"\end{figure}"+"\n")
                f.write(r"\clearpage"+"\n")
            f.write(r"\end{document}")
        os.system('pdflatex -interaction=batchmode {0}'.format(fname))
        self.logger.info("Done.")

    def lineplot(self):
        """
        Generate PDF of channel cube diagnostic plots

        Inputs: Nothing

        Returns: Nothing
        """
        #
        # Setup parameters for zooming
        #
        center_x = int(self.cp['imsize'][0]/2)
        center_y = int(self.cp['imsize'][1]/2)
        #
        # Get center pixels
        #
        goodplots = []
        for spw in self.line_spws.split(','):
            # check that this spectral window exists
            fname = '{0}.spw{1}.channel.clean.image.fits'.format(self.field, spw)
            if self.uvtaper:
                fname = '{0}.spw{1}.channel.clean.uvtaper.image.fits'.format(self.field, spw)
            if not os.path.exists(fname):
                continue
            #
            # Loop over all plot filenames
            #
            if self.uvtaper:
                fitsfiles = ['{0}.spw{1}.channel.dirty.uvtaper.image.fits'.format(self.field, spw),
                             '{0}.spw{1}.channel.clean.uvtaper.image.fits'.format(self.field, spw),
                             '{0}.spw{1}.channel.clean.uvtaper.residual.fits'.format(self.field, spw),
                             '{0}.spw{1}.channel.clean.uvtaper.pbcor.image.fits'.format(self.field, spw)]
                maskfile = '{0}.spw{1}.channel.clean.uvtaper.mask.fits'.format(self.field, spw)
                titles = ['{0} - spw: {1} - Dirty/UVTap'.format(self.field, spw),
                          '{0} - spw: {1} - Clean/UVTap'.format(self.field, spw),
                          '{0} - spw: {1} - Residual/UVTap'.format(self.field, spw),
                          '{0} - spw: {1} - UVTap/PBCorr'.format(self.field, spw)]
                labels = ['Flux Density (Jy/beam)',
                          'Flux Density (Jy/beam)',
                          'Flux Density (Jy/beam)',
                          'Flux Density (Jy/beam)']
                vlims = [(None, None),
                         (None, None),
                         (None, None),
                         (None, None)]
            else:
                fitsfiles = ['{0}.spw{1}.channel.dirty.image.fits'.format(self.field, spw),
                             '{0}.spw{1}.channel.clean.image.fits'.format(self.field, spw),
                             '{0}.spw{1}.channel.clean.residual.fits'.format(self.field, spw),
                             '{0}.spw{1}.channel.clean.pbcor.image.fits'.format(self.field, spw)]
                maskfile = '{0}.spw{1}.channel.clean.mask.fits'.format(self.field, spw)
                titles = ['{0} - spw: {1} - Dirty'.format(self.field, spw),
                          '{0} - spw: {1} - Clean'.format(self.field, spw),
                          '{0} - spw: {1} - Residual'.format(self.field, spw),
                          '{0} - spw: {1} - PBCorr'.format(self.field, spw)]
                labels = ['Flux Density (Jy/beam)',
                          'Flux Density (Jy/beam)',
                          'Flux Density (Jy/beam)',
                          'Flux Density (Jy/beam)']
                vlims = [(None, None),
                         (None, None),
                         (None, None),
                         (None, None)]
            for fitsfile, title, label, vlim in \
                zip(fitsfiles, titles, labels, vlims):
                #
                # Open fits file, generate WCS
                #
                hdulist = fits.open(fitsfile)
                hdu = hdulist[0]
                wcs = WCS(hdu.header)
                #
                # Generate figure
                #
                plt.ioff()
                fig = plt.figure()
                ax = plt.subplot(projection=wcs.sub(['celestial']))
                ax.set_title(title)
                center_chan = hdu.data.shape[1]/2
                img = ax.imshow(hdu.data[0,center_chan],
                                origin='lower', interpolation='none',
                                cmap='binary', vmin=vlim[0],
                                vmax=vlim[1])
                # ax.grid(True,color='black',ls='solid')
                ax.coords[0].set_major_formatter('hh:mm:ss')
                ax.set_xlabel('RA (J2000)')
                ax.set_ylabel('Declination (J2000)')
                #
                # Plot clean mask as contour
                #
                maskhdu = fits.open(maskfile)[0]
                contour = ax.contour(maskhdu.data[0, center_chan],
                                     levels=[0.5], origin='lower',
                                     colors='r', linewidths=2)
                #
                # Plot beam, if it is defined
                #
                if 'BMAJ' in hdu.header.keys():
                    cell = float(self.cp['cell'].replace('arcsec', ''))
                    beam_maj = hdu.header['BMAJ']*3600./cell
                    beam_min = hdu.header['BMIN']*3600./cell
                    beam_pa = hdu.header['BPA']
                    ellipse = Ellipse((center_x-int(3.*center_x/4),
                                       center_y-int(3.*center_y/4)),
                                      beam_min, beam_maj,
                                      angle=beam_pa, fill=True,
                                      zorder=10, hatch='///',
                                      edgecolor='black',
                                      facecolor='white')
                    ax.add_patch(ellipse)
                elif len(hdulist) > 1:
                    hdu = hdulist[1]
                    cell = float(self.cp['cell'].replace('arcsec',''))
                    beam_maj = hdu.data['BMAJ'][center_chan]/cell
                    beam_min = hdu.data['BMIN'][center_chan]/cell
                    beam_pa = hdu.data['BPA'][center_chan]
                    ellipse = Ellipse((center_x-int(3.*center_x/4),
                                       center_y-int(3.*center_y/4)),
                                      beam_min, beam_maj,
                                      angle=beam_pa, fill=True,
                                      zorder=10, hatch='///',
                                      edgecolor='black',
                                      facecolor='white')
                    ax.add_patch(ellipse)
                #
                # Plot colorbar
                #
                cbar = fig.colorbar(img, fraction=0.046, pad=0.04)
                cbar.set_label(label)
                #
                # Re-scale to fit, then save
                #
                fig.savefig(fitsfile.replace('.fits', '.pdf'),
                            bbox_inches='tight')
                plt.close(fig)
                plt.ion()
                goodplots.append(fitsfile.replace('.fits', '.pdf'))
            #
            # Generate spectrum
            #
            if self.uvtaper:
                fitsfile = '{0}.spw{1}.channel.clean.uvtaper.image.fits'.format(self.field, spw)
            else:
                fitsfile = '{0}.spw{1}.channel.clean.image.fits'.format(self.field, spw)
            hdu = fits.open(fitsfile)[0]
            spec = hdu.data[0, :, center_x, center_y]
            isnan = spec == 0.
            spec[isnan] = np.nan
            velo = (np.arange(len(spec))*hdu.header['CDELT3'] + hdu.header['CRVAL3'])/1000.
            #
            # Generate figure
            #
            plt.ioff()
            fig = plt.figure()
            ax = plt.subplot()
            ax.plot(velo, spec, 'k-')
            ax.set_xlabel('Velocity (km/s)')
            ax.set_ylabel('Flux Density (Jy/beam)')
            ax.set_xlim(np.nanmin(velo), np.nanmax(velo))
            ybuff = 0.1*(np.nanmax(spec)-np.nanmin(spec))
            ax.set_ylim(np.nanmin(spec)-ybuff, np.nanmax(spec)+ybuff)
            if self.uvtaper:
                ax.set_title('{0} - spw {1} - UVTap/Center'.format(self.field, spw))
            else:
                ax.set_title('{0} - spw {1} - Center'.format(self.field, spw))
            ax.grid(False)
            fig.tight_layout()
            fig.savefig(fitsfile.replace('.fits', '.spec.pdf'),
                        bbox_inches='tight')
            plt.close(fig)
            plt.ion()
            goodplots.append(fitsfile.replace('.fits', '.spec.pdf'))
        #
        # Generate PDF of plots
        #
        # need to fix filenames so LaTeX doesn't complain
        goodplots = ['{'+fn.replace('.pdf', '')+'}.pdf' for fn in goodplots]
        self.logger.info("Generating PDF...")
        fname = '{0}.lineplots.tex'.format(self.field)
        if self.uvtaper:
            fname = '{0}.uvtaper.lineplots.tex'.format(self.field)
        with open(fname, 'w') as f:
            f.write(r"\documentclass{article}"+"\n")
            f.write(r"\usepackage{graphicx}"+"\n")
            f.write(r"\usepackage[margin=0.1cm]{geometry}"+"\n")
            f.write(r"\begin{document}"+"\n")
            for i in range(0, len(goodplots), 5):
                f.write(r"\begin{figure}"+"\n")
                f.write(r"\centering"+"\n")
                f.write(r"\includegraphics[width=0.45\textwidth]{"+goodplots[i]+"}\n")
                f.write(r"\includegraphics[width=0.45\textwidth]{"+goodplots[i+1]+r"} \\"+"\n")
                f.write(r"\includegraphics[width=0.45\textwidth]{"+goodplots[i+2]+"}\n")
                f.write(r"\includegraphics[width=0.45\textwidth]{"+goodplots[i+3]+r"} \\"+"\n")
                f.write(r"\includegraphics[width=0.45\textwidth]{"+goodplots[i+4]+"}\n")
                f.write(r"\end{figure}"+"\n")
                f.write(r"\clearpage"+"\n")
            f.write(r"\end{document}")
        os.system('pdflatex -interaction=batchmode {0}'.format(fname))
        self.logger.info("Done.")

def main(vis, field, config_file, stokes='I', spws='',
         uvtaper=False, interactive=False, savemodel=False, auto=''):
    """
    Generate and clean images

    Inputs:
      vis :: string
        The measurement set containing all data for field
      field :: string
        The field name to image
      config_file :: string
        filename of the configuration file for this project
      stokes :: string
        The Stokes parameters we're imaging. e.g. 'I' or 'IQUV'
      spws :: string
        comma-separated list of spws to clean
        if empty, clean all spws
      uvtaper :: boolean
        if True, apply UV tapering
      interactive :: boolean
        if True, interactively clean
      savemodel :: boolean
        if True, save individual MFS images of each spectral window
        to the model column of the measurement set for
        self-calibration. This can only be done with stokes='I'
      auto :: string
        if not an empty string, it is a comma separated
        list of menu items to perform, i.e. auto='0,1,4,5,6'

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # Check inputs
    #
    if not os.path.isdir(vis):
        logger.critical('Measurement set not found!')
        raise ValueError('Measurement set not found!')
    if not os.path.exists(config_file):
        logger.critical('Configuration file not found')
        raise ValueError('Configuration file not found!')
    if savemodel and stokes != 'I':
        logger.critical('Can only save visibility model with Stokes I')
        raise ValueError('Can only save visibility model with Stokes I')
    #
    # load configuration file
    #
    config = ConfigParser.ConfigParser()
    logger.info("Reading configuration file {0}".format(config_file))
    config.read(config_file)
    logger.info("Done.")
    #
    # Initialize Imaging object
    #
    imag = Imaging(vis, field, logger, config, uvtaper=uvtaper,
                   spws=spws, stokes=stokes, savemodel=savemodel,
                   interactive=interactive)
    #
    # Prompt the user with a menu for each option, or auto-do them
    #
    auto_items = auto.split(',')
    auto_ind = 0
    while True:
        if not auto:
            print("0. Dirty image combined continuum spws (MFS; multi-term; multi-scale)")
            print("1. Clean combined continuum spws (MFS; multi-term; multi-scale)")
            print("2. Dirty image each continuum spw (MFS; multi-scale)")
            print("3. Clean each continuum spw (MFS; multi-scale)")
            print("4. Dirty image each continuum spw (channel; multi-scale)")
            print("5. Clean each continuum spw (channel; multi-scale)")
            print("6. Dirty image each line spw (MFS; multi-scale)")
            print("7. Clean each line spw (MFS; multi-scale)")
            print("8. Dirty image each line spw (channel; multi-scale)")
            print("9. Clean each line spw (channel; multi-scale)")
            print("10. Generate continuum diagnostic plots")
            print("11. Generate spectral line diagnostic plots")
            print("q [quit]")
            answer = input("> ")
        else:
            answer = auto_items[auto_ind]
            auto_ind += 1
        if answer == '0':
            imag.mfs_dirty_cont()
        elif answer == '1':
            imag.mfs_clean_cont()
        elif answer == '2':
            imag.mfs_dirty_spws(imag.cont_spws)
        elif answer == '3':
            imag.mfs_clean_spws(imag.cont_spws)
        elif answer == '4':
            imag.channel_dirty_spws(imag.cont_spws, 'cont')
        elif answer == '5':
            imag.channel_clean_spws(imag.cont_spws, 'cont')
        elif answer == '6':
            imag.mfs_dirty_spws(imag.line_spws)
        elif answer == '7':
            imag.mfs_clean_spws(imag.line_spws)
        elif answer == '8':
            imag.channel_dirty_spws(imag.line_spws, 'line')
        elif answer == '9':
            imag.channel_clean_spws(imag.line_spws, 'line')
        elif answer == '10':
            imag.contplot()
        elif answer == '11':
            imag.lineplot()
        elif answer.lower() == 'q' or answer.lower() == 'quit':
            break
        else:
            print("Input not recognized.")
        if auto_ind >= len(auto_items):
            break
