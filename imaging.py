"""
imaging.py - WISP Imaging Pipeline

Imaging a single field in a measurement set by scripting CASA tasks
and generating diagnostic plots.

Copyright(C) 2018 by
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
"""

import __main__ as casa
import os
import gc
import numpy as np
import glob
import logging
import logging.config
import ConfigParser
import shutil
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.patches import Ellipse

__VERSION__ = "1.0"

# load logging configuration file
logging.config.fileConfig('logging.conf')

def natural_sort(l):
    """
    Natural sort an alphanumeric list

    Inputs: l
      l :: list of strings
        The list of strings to be sorted

    Returns: sorted_l
      sorted_l :: list of strings
        The sorted list
    """
    # Convert text to integers or lowercase
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    # define the sorting algorithm
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    # return the sorted list
    return sorted(l, key = alphanum_key)

def setup(vis='', uvtaper=False, config=None):
    """
    Perform setup tasks: find line and continuum spectral windows
                         get clean parameters

    Inputs:
      vis :: string
        The measurement set
      uvtaper :: boolean
        If True, get UV tapering clean parameters
      config :: a ConfigParser object
        The ConfigParser object for this project

    Returns: my_cont_spws, my_line_spws, cp
      my_cont_spws :: string
        comma-separated string of continuum spws
      my_line_spws :: string
        comma-separated string of line spws
      cp :: dictionary
        clean parameters read from config file
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # check config
    #
    if config is None:
        logger.critical("Error: Need to supply a config")
        raise ValueError("Config is None") 
    #
    # Get continuum and line spws from configuration file
    #
    my_cont_spws = config.get("Spectral Windows", "Continuum")
    my_line_spws = config.get("Spectral Windows", "Line")
    logger.info("Found continuum spws: {0}".format(my_cont_spws))
    logger.info("Found line spws: {0}".format(my_line_spws))
    #
    # Get clean parameters from configuration file
    #
    cp = {}
    cp["lineids"] = config.get("Clean", "lineids").split(',')
    cp["restfreqs"] = config.get("Clean", "restfreqs").split(',')
    # clean info
    cp["imsize"] = [int(foo) for foo in config.get("Clean","imsize").split(',')]
    cp["pblimit"] = config.getfloat("Clean", "pblimit")
    cp["cell"] = "{0}arcsec".format(config.getfloat("Clean", "cell"))
    cp["weighting"] = config.get("Clean", "weighting")
    cp["robust"] = config.getfloat("Clean", "robust")
    cp["scales"] = [int(foo) for foo in config.get("Clean", "scales").split(',') if foo != '']
    cp["gain"] = config.getfloat("Clean", "gain")
    cp["cyclefactor"] = config.getfloat("Clean", "cyclefactor")
    cp["lightniter"] = config.getint("Clean", "lightniter")
    cp["maxniter"] = config.getint("Clean", "maxniter")
    cp["nrms"] = config.getfloat("Clean", "nrms")
    cp["contpbchan"] = config.get("Clean", "contpbchan")
    cp["nterms"] = config.getint("Clean", "nterms")
    if uvtaper:
        cp["outertaper"] = ["{0}arcsec".format(config.getfloat("Clean","outertaper"))]
    else:
        cp["outertaper"] = []
    # re-grid info
    velstart = config.getfloat("Clean", "start")
    cp["velstart"] = "{0}km/s".format(velstart)
    chanwidth = config.getfloat("Clean", "width")
    cp["chanwidth"] = "{0}km/s".format(chanwidth)
    cp["nchan"] = config.getint("Clean", "nchan")
    cp["chanbuffer"] = config.getint("Clean", "chanbuffer")
    cp["cvelstart"] = "{0}km/s".format(velstart-(cp["chanbuffer"]*chanwidth))
    cp["cvelnchan"] = cp["nchan"]+2*cp["chanbuffer"]
    cp["lineoutframe"] = config.get("Clean", "lineoutframe")
    cp["veltype"] = config.get("Clean", "veltype")
    cp["interpolation"] = config.get("Clean", "interpolation")
    cp["contoutframe"] = config.get("Clean", "contoutframe")
    if uvtaper:
        heading = "Mask Taper"
    else:
        heading = "Mask NoTaper"
    cp["contpbmask"] = config.getfloat(heading, "contpbmask")
    cp["contsidelobethreshold"] = config.getfloat(heading, "contsidelobethreshold")
    cp["contnoisethreshold"] = config.getfloat(heading, "contnoisethreshold")
    cp["contlownoisethreshold"] = config.getfloat(heading, "contlownoisethreshold")
    cp["contnegativethreshold"] = config.getfloat(heading, "contnegativethreshold")
    cp["contsmoothfactor"] = config.getfloat(heading, "contsmoothfactor")
    cp["contminbeamfrac"] = config.getfloat(heading, "contminbeamfrac")
    cp["contcutthreshold"] = config.getfloat(heading, "contcutthreshold")
    cp["contgrowiterations"] = config.getint(heading, "contgrowiterations")
    cp["linepbmask"] = config.getfloat(heading, "linepbmask")
    cp["linesidelobethreshold"] = config.getfloat(heading, "linesidelobethreshold")
    cp["linenoisethreshold"] = config.getfloat(heading, "linenoisethreshold")
    cp["linelownoisethreshold"] = config.getfloat(heading, "linelownoisethreshold")
    cp["linenegativethreshold"] = config.getfloat(heading, "linenegativethreshold")
    cp["linesmoothfactor"] = config.getfloat(heading, "linesmoothfactor")
    cp["lineminbeamfrac"] = config.getfloat(heading, "lineminbeamfrac")
    cp["linecutthreshold"] = config.getfloat(heading, "linecutthreshold")
    cp["linegrowiterations"] = config.getint(heading, "linegrowiterations")
    return (my_cont_spws, my_line_spws, cp)

def regrid_velocity(vis='', spws='', cp={}, config=None):
    """
    Re-grid velocity axis of each spectral window

    Inputs:
      vis :: string
        The measurement set
      spws :: string
        comma-separated list of spectral windows to regrid
      cp :: dictionary
        clean parameters
      config :: a ConfigParser object
        The ConfigParser object for this project

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # check config
    #
    if config is None:
        logger.critical("Error: Need to supply a config")
        raise ValueError("Config is None") 
    #
    # Re-grid velocity axis of line spectral windows
    #
    for spw in spws.split(','):
        regrid_vis = vis.replace('.ms','')+'.spw{0}.cvel.ms'.format(spw)
        #
        # Check if it already exists
        #
        if os.path.isdir(regrid_vis):
            logger.info("Found {0}".format(regrid_vis))
            continue
        #
        # If not, re-grid it
        #
        logger.info("Regridding velocity axis of spw {0}".format(spw))
        spw_ind = config.get("Spectral Windows","Line").split(',').index(spw)
        restfreq = config.get("Clean","restfreqs").split(',')[spw_ind]
        casa.cvel2(vis=vis, outputvis=regrid_vis, spw=spw,
                   restfreq=restfreq, mode='velocity', 
                   start=cp['cvelstart'], width=cp['chanwidth'], 
                   nchan=cp['cvelnchan'], outframe=cp['lineoutframe'], 
                   veltype=cp['veltype'], interpolation=cp['interpolation'])
        logger.info("Done.")

def mfs_dirty_cont(field='', vis='', my_cont_spws='', cp={},
                   uvtaper=False):
    """
    Dirty image continuum spws using multi-frequency synthesis

    Inputs:
      field :: string
        The field to be cleaned
      vis :: string
        The measurement set
      my_cont_spws :: string
        comma-separated string of continuum spws
      cp :: dictionary
        dictionary of clean parameters
      uvtaper :: boolean
        if True, apply UV tapering

    Returns:
      Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # Dirty image continuum
    #
    imagename='{0}.cont.mfs.dirty'.format(field)
    if uvtaper:
        imagename = imagename + '.uvtaper'
    logger.info("Generating dirty continuum image (MFS)...")
    casa.tclean(vis=vis, imagename=imagename, field=field, spw=my_cont_spws, 
                specmode='mfs', threshold='0mJy', niter=0, 
                nterms=cp['nterms'], deconvolver='mtmfs', scales=cp['scales'], 
                gain=cp['gain'], cyclefactor=cp['cyclefactor'], 
                imsize=cp['imsize'], pblimit=-1.0, cell=cp['cell'], 
                weighting=cp['weighting'], robust=cp['robust'], 
                uvtaper=cp['outertaper'], pbcor=False)
    logger.info("Done.")
    #
    # Primary beam correction using PB of center channel
    #
    logger.info("Performing primary beam correction...")
    spwlist = [int(spw) for spw in my_cont_spws.split(',')]
    weightlist = [1.0 for spw in spwlist]
    chanlist = [int(cp['contpbchan']) for foo in my_cont_spws.split(',')]
    casa.widebandpbcor(vis=vis, imagename=imagename, 
                       nterms=cp['nterms'], pbmin=cp['pblimit'], threshold='0.1mJy', 
                       spwlist=spwlist, weightlist=weightlist, chanlist=chanlist)
    logger.info("Done.")
    #
    # Export to fits
    #
    logger.info("Exporting fits file...")
    casa.exportfits(imagename='{0}.pbcor.workdirectory/{0}.pb.tt0'.format(imagename), 
                    fitsimage='{0}.cont.mfs.pb.fits'.format(field), 
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
    logger.info("Done.")

def mfs_clean_cont(field='', vis='', my_cont_spws='', cp={}, 
                   uvtaper=False, interactive=False):
    """
    Clean continuum spws using multi-frequency synthesis

    Inputs:
      field :: string
        The field to be cleaned
      vis :: string
        The measurement set
      my_cont_spws :: string
        comma-separated string of continuum spws
      cp :: dictionary
        dictionary of clean parameters
      uvtaper :: boolean
        if True, apply UV tapering
      interactive :: boolean
        if True, clean interactively

    Returns:
      Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # Lightly clean continuum to get RMS threshold
    #
    imagename='{0}.cont.mfs.lightclean'.format(field)
    if uvtaper:
        imagename = imagename + '.uvtaper'
    logger.info("Lightly cleaning continuum image (MFS)...")
    casa.tclean(vis=vis, imagename=imagename, field=field, spw=my_cont_spws, 
                specmode='mfs', threshold='0mJy', niter=cp["lightniter"], 
                usemask='auto-multithresh', pbmask=cp['contpbmask'], 
                sidelobethreshold=cp['contsidelobethreshold'], 
                noisethreshold=cp['contnoisethreshold'], 
                lownoisethreshold=cp['contlownoisethreshold'], 
                negativethreshold=cp['contnegativethreshold'], 
                smoothfactor=cp['contsmoothfactor'], minbeamfrac=cp['contminbeamfrac'], 
                cutthreshold=cp['contcutthreshold'], growiterations=cp['contgrowiterations'], 
                nterms=cp['nterms'], deconvolver='mtmfs', 
                scales=cp['scales'], 
                gain=cp['gain'], cyclefactor=cp['cyclefactor'], 
                imsize=cp['imsize'], pblimit=-1.0, cell=cp['cell'], 
                weighting=cp['weighting'], robust=cp['robust'], 
                uvtaper=cp['outertaper'], pbcor=False)
    logger.info("Done.")
    #
    # Get RMS of residuals outside of clean mask
    #
    dat = casa.imstat(imagename='{0}.residual.tt0'.format(imagename), 
                      axes=[0, 1], mask="'{0}.mask' == 0".format(imagename))
    logger.info("Max un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.max(dat['rms'])))
    logger.info("Mean un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.mean(dat['rms'])))
    logger.info("Median un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.median(dat['rms'])))
    logger.info("Max un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.max(dat['medabsdevmed'])))
    logger.info("Mean un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.mean(dat['medabsdevmed'])))
    logger.info("Median un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.median(dat['medabsdevmed'])))
    logger.info("Using median MADS*1.4826 times {0} (user-defined) as threshold".format(cp['nrms']))
    threshold = '{0:.2f}mJy'.format(cp["nrms"]*1000.*1.4826*np.median(dat['medabsdevmed']))
    #
    # Clean to threshold
    #
    imagename='{0}.cont.mfs.clean'.format(field)
    if uvtaper:
        imagename = imagename + '.uvtaper'
    logger.info("Cleaning continuum image (MFS) to threshold: {0}...".format(threshold))
    casa.tclean(vis=vis, imagename=imagename, field=field, spw=my_cont_spws, 
                specmode='mfs', threshold=threshold, niter=cp['maxniter'], 
                usemask='auto-multithresh', pbmask=cp['contpbmask'], 
                sidelobethreshold=cp['contsidelobethreshold'], 
                noisethreshold=cp['contnoisethreshold'], 
                lownoisethreshold=cp['contlownoisethreshold'], 
                negativethreshold=cp['contnegativethreshold'], 
                smoothfactor=cp['contsmoothfactor'], minbeamfrac=cp['contminbeamfrac'], 
                cutthreshold=cp['contcutthreshold'], growiterations=cp['contgrowiterations'], 
                nterms=cp['nterms'], deconvolver='mtmfs', 
                scales=cp['scales'], 
                gain=cp['gain'], cyclefactor=cp['cyclefactor'], 
                imsize=cp['imsize'], pblimit=-1.0, cell=cp['cell'], 
                weighting=cp['weighting'], robust=cp['robust'], 
                uvtaper=cp['outertaper'], pbcor=False, 
                interactive=interactive)
    logger.info("Done.")
    #
    # Primary beam correction using PB of center channel
    #
    logger.info("Performing primary beam correction...")
    spwlist = [int(spw) for spw in my_cont_spws.split(',')]
    weightlist = [1.0 for spw in spwlist]
    chanlist = [int(cp['contpbchan']) for foo in my_cont_spws.split(',')]
    casa.widebandpbcor(vis=vis, imagename=imagename, 
                       nterms=cp['nterms'], pbmin=cp['pblimit'], threshold='0.1mJy', 
                       spwlist=spwlist, weightlist=weightlist, chanlist=chanlist)
    logger.info("Done.")
    #
    # Export to fits
    #
    logger.info("Exporting fits file...")
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
    logger.info("Done.")

def mfs_dirty_spws(field='', vis='', my_spws='',
                   cp={}, uvtaper=False):
    """
    Dirty image each supplied spw (MFS)

    Inputs:
      field :: string
        The field to be imaged
      vis :: string
        The measurement set
      my_spws :: string
        comma-separated string of spws to image
      cp :: dictionary
        dictionary of clean parameters
      uvtaper :: boolean
        if True, apply UV tapering

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    for spw in my_spws.split(','):
        #
        # Dirty image
        #
        imagename='{0}.spw{1}.mfs.dirty'.format(field,spw)
        if uvtaper:
            imagename = imagename + '.uvtaper'
        logger.info("Generating dirty image of spw {0} (MFS)...".format(spw))
        casa.tclean(vis=vis, imagename=imagename, field=field, spw=spw, 
                    specmode='mfs', threshold='0mJy', niter=0, 
                    deconvolver='multiscale', scales=cp['scales'], 
                    gain=cp['gain'], cyclefactor=cp['cyclefactor'], 
                    imsize=cp['imsize'], pblimit=-1.0, cell=cp['cell'], 
                    weighting=cp['weighting'], robust=cp['robust'], 
                    uvtaper=cp['outertaper'], pbcor=False)
        logger.info("Done.")
        #
        # Generate primary beam image
        #
        pbimagename='{0}.spw{1}.mfs.pb'.format(field, spw)
        if uvtaper:
            pbimagename = pbimagename + '.uvtaper'
        logger.info("Generating primary beam image of spw {0} (MFS)...".format(spw))
        casa.tclean(vis=vis, imagename=pbimagename, field=field, spw=spw, 
                    specmode='mfs', threshold='0mJy', niter=0, 
                    deconvolver='multiscale', scales=cp['scales'], 
                    gain=cp['gain'], cyclefactor=cp['cyclefactor'], 
                    imsize=cp['imsize'], pblimit=cp['pblimit'], cell=cp['cell'], 
                    weighting=cp['weighting'], robust=cp['robust'], 
                    uvtaper=cp['outertaper'], pbcor=False)
        logger.info("Done.")
        #
        # Primary beam correction
        #
        pbimagename='{0}.spw{1}.mfs.pb'.format(field, spw)
        if uvtaper:
            pbimagename = pbimagename + '.uvtaper'
        logger.info("Performing primary beam correction...")
        casa.impbcor(imagename='{0}.image'.format(imagename), 
                     pbimage='{0}.pb'.format(pbimagename), 
                     outfile='{0}.pbcor.image'.format(imagename))
        logger.info("Done.")
        #
        # Export to fits
        #
        logger.info("Exporting fits file...")
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
        logger.info("Done.")

def mfs_clean_spws(field='', vis='', my_spws='',
                   cp={}, uvtaper=False, interactive=False):
    """
    Clean each supplied spw (MFS)

    Inputs:
      field :: string
        The field to be imaged
      vis :: string
        The measurement set
      my_spws :: string
        comma-separated string of spws to image
      cp :: dictionary
        dictionary of clean parameters
      uvtaper :: boolean
        if True, apply UV tapering
      interactive :: boolean
        if True, clean interactively

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    for spw in my_spws.split(','):
        #
        # Lightly clean to get threshold
        #
        imagename='{0}.spw{1}.mfs.lightclean'.format(field,spw)
        if uvtaper:
            imagename = imagename + '.uvtaper'
        logger.info("Lightly cleaning spw {0} (MFS)...".format(spw))
        casa.tclean(vis=vis, imagename=imagename, field=field, spw=spw, 
                    specmode='mfs', threshold='0mJy', niter=cp["lightniter"], 
                    usemask='auto-multithresh', pbmask=cp['linepbmask'], 
                    sidelobethreshold=cp['linesidelobethreshold'], 
                    noisethreshold=cp['linenoisethreshold'], 
                    lownoisethreshold=cp['linelownoisethreshold'], 
                    negativethreshold=cp['linenegativethreshold'], 
                    smoothfactor=cp['linesmoothfactor'], minbeamfrac=cp['lineminbeamfrac'], 
                    cutthreshold=cp['linecutthreshold'], growiterations=cp['linegrowiterations'], 
                    deconvolver='multiscale', scales=cp['scales'], 
                    gain=cp['gain'], cyclefactor=cp['cyclefactor'], 
                    imsize=cp['imsize'], pblimit=-1.0, cell=cp['cell'], 
                    weighting=cp['weighting'], robust=cp['robust'], 
                    uvtaper=cp['outertaper'], pbcor=False)
        logger.info("Done.")
        #
        # Get RMS of residuals
        #
        dat = casa.imstat(imagename='{0}.residual'.format(imagename), 
                          axes=[0, 1], mask="'{0}.mask' == 0".format(imagename))
        logger.info("Max un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.max(dat['rms'])))
        logger.info("Mean un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.mean(dat['rms'])))
        logger.info("Median un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.median(dat['rms'])))
        logger.info("Max un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.max(dat['medabsdevmed'])))
        logger.info("Mean un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.mean(dat['medabsdevmed'])))
        logger.info("Median un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.median(dat['medabsdevmed'])))
        logger.info("Using median MADS*1.4826 times {0} (user-defined) as threshold".format(cp['nrms']))
        threshold = '{0:.2f}mJy'.format(cp["nrms"]*1000.*1.4826*np.median(dat['medabsdevmed']))
        #
        # Clean to threshold
        #
        imagename='{0}.spw{1}.mfs.clean'.format(field, spw)
        if uvtaper:
            imagename = imagename + '.uvtaper'
        logger.info("Cleaning spw {0} (MFS) to threshold: {1}...".format(spw, threshold))
        casa.tclean(vis=vis, imagename=imagename, field=field, spw=spw, 
                    specmode='mfs', threshold=threshold, niter=cp["maxniter"], 
                    usemask='auto-multithresh', pbmask=cp['linepbmask'], 
                    sidelobethreshold=cp['linesidelobethreshold'], 
                    noisethreshold=cp['linenoisethreshold'], 
                    lownoisethreshold=cp['linelownoisethreshold'], 
                    negativethreshold=cp['linenegativethreshold'], 
                    smoothfactor=cp['linesmoothfactor'], minbeamfrac=cp['lineminbeamfrac'], 
                    cutthreshold=cp['linecutthreshold'], growiterations=cp['linegrowiterations'], 
                    deconvolver='multiscale', scales=cp['scales'], 
                    gain=cp['gain'], cyclefactor=cp['cyclefactor'], 
                    imsize=cp['imsize'], pblimit=-1.0, cell=cp['cell'], 
                    weighting=cp['weighting'], robust=cp['robust'], 
                    uvtaper=cp['outertaper'], pbcor=False, 
                    interactive=interactive)
        logger.info("Done.")
        #
        # Primary beam correction
        #
        pbimagename='{0}.spw{1}.mfs.pb'.format(field, spw)
        if uvtaper:
            pbimagename = pbimagename + '.uvtaper'
        logger.info("Performing primary beam correction...")
        casa.impbcor(imagename='{0}.image'.format(imagename), 
                     pbimage='{0}.pb'.format(pbimagename), 
                     outfile='{0}.pbcor.image'.format(imagename))
        logger.info("Done.")
        #
        # Export to fits
        #
        logger.info("Exporting fits file...")
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
        logger.info("Done.")

def channel_dirty_spws(field='', vis='', my_spws='',
                       cp={}, spwtype='line', regrid=False,
                       uvtaper=False, config=None):
    """
    Dirty image each supplied spectral window (channel cube)

    Inputs:
      field :: string
        The field to be imaged
      vis :: string
        The measurement set
      my_spws :: string
        comma-separated string of spws to image
      cp :: dictionary
        dictionary of clean parameters
      spwtype :: string
        'line' or 'cont', determine which clean parameters to use
      regrid :: boolean
        if True, look for CVEL2 measurement set
      uvtaper :: boolean
        if True, apply UV tapering
      config :: a ConfigParser object
        ConfigParser object for this project

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # check config
    #
    if config is None:
        logger.critical("Error: Need to supply a config")
        raise ValueError("Config is None") 
    #
    # Set channel parameters
    #
    if spwtype == 'cont':
        restfreqs = [None for spw in my_spws.split(',')]
        start = None
        width = None
        nchan = None
        outframe = cp['contoutframe']
        veltype = None
        interpolation = None
    elif spwtype == 'line':
        spw_inds = [config.get("Spectral Windows","Line").split(',').index(spw) for spw in my_spws.split(',')]
        restfreqs = [config.get("Clean","restfreqs").split(',')[spw_ind] for spw_ind in spw_inds]
        start = cp['velstart']
        width = cp['chanwidth']
        nchan = cp['nchan']
        outframe = cp['lineoutframe']
        veltype = cp['veltype']
        interpolation = cp['interpolation']
    else:
        logger.critical("Error: spwtype {0} not supported".format(spwtype))
        raise ValueError("Invalid spwtype")
    #
    # Loop over spws
    #
    for spw,restfreq in zip(my_spws.split(','),restfreqs):
        #
        # dirty image spw
        #
        if regrid:
            myvis = vis.replace('.ms','')+'.spw{0}.cvel.ms'.format(spw)
            myspw = '0'
        else:
            myvis = vis
            myspw = spw
        imagename='{0}.spw{1}.channel.dirty'.format(field,spw)
        if uvtaper:
            imagename = imagename + '.uvtaper'
        logger.info("Dirty imaging spw {0} (restfreq: {1})...".format(spw,restfreq))
        casa.tclean(vis=myvis, imagename=imagename, field=field, spw=myspw, 
                    specmode='cube', threshold='0mJy', niter=0, 
                    deconvolver='multiscale', scales=cp['scales'], 
                    gain=cp['gain'], cyclefactor=cp['cyclefactor'], 
                    imsize=cp['imsize'], pblimit=-1.0, cell=cp['cell'], 
                    weighting=cp['weighting'], robust=cp['robust'], 
                    restfreq=restfreq, start=start, width=width, nchan=nchan, 
                    outframe=outframe, veltype=veltype, interpolation=interpolation, 
                    uvtaper=cp['outertaper'], pbcor=False)
        logger.info("Done.")
        #
        # Generate primary beam image
        #
        pbimagename='{0}.spw{1}.channel.pb'.format(field, spw)
        if uvtaper:
            pbimagename = pbimagename + '.uvtaper'
        logger.info("Generating primary beam image of spw {0} (MFS)...".format(spw))
        casa.tclean(vis=myvis, imagename=pbimagename, field=field, spw=myspw, 
                    specmode='cube', threshold='0mJy', niter=0, 
                    deconvolver='multiscale', scales=cp['scales'], 
                    gain=cp['gain'], cyclefactor=cp['cyclefactor'], 
                    imsize=cp['imsize'], pblimit=cp['pblimit'], cell=cp['cell'], 
                    weighting=cp['weighting'], robust=cp['robust'], 
                    restfreq=restfreq, start=start, width=width, nchan=nchan, 
                    outframe=outframe, veltype=veltype, interpolation=interpolation, 
                    uvtaper=cp['outertaper'], pbcor=False)
        logger.info("Done.")
        #
        # Primary beam correction
        #
        pbimagename='{0}.spw{1}.channel.pb'.format(field, spw)
        if uvtaper:
            pbimagename = pbimagename + '.uvtaper'
        logger.info("Performing primary beam correction...")
        casa.impbcor(imagename='{0}.image'.format(imagename), 
                     pbimage='{0}.pb'.format(pbimagename), 
                     outfile='{0}.pbcor.image'.format(imagename))
        logger.info("Done.")
        #
        # Export to fits
        #
        logger.info("Exporting fits file...")
        casa.exportfits(imagename='{0}.image'.format(imagename), 
                        fitsimage='{0}.image.fits'.format(imagename), 
                        velocity=True, overwrite=True, history=False)
        casa.exportfits(imagename='{0}.residual'.format(imagename), 
                        fitsimage='{0}.residual.fits'.format(imagename), 
                        velocity=True, overwrite=True, history=False)
        casa.exportfits(imagename='{0}.pbcor.image'.format(imagename), 
                        fitsimage='{0}.pbcor.image.fits'.format(imagename), 
                        velocity=True, overwrite=True, history=False)
        logger.info("Done.")

def channel_clean_spws(field='', vis='', my_spws='', 
                       cp={}, spwtype='line', regrid=False, 
                       uvtaper=False, interactive=False, config=None):
    """
    Clean all supplied spws by channel using clean mask from MFS
    images.

    Inputs:
      field :: string
        The field to be imaged
      vis :: string
        The measurement set
      my_spws :: string
        comma-separated string of spws to image
      cp :: dictionary
        dictionary of clean parameters
      spwtype :: string
        'line' or 'cont', determine which clean parameters to use
      regrid :: boolean
         if True, look for CVEL2 measurement set
      uvtaper :: boolean
        if True, apply UV tapering
      interactive :: boolean
        if True, clean interactively
      config :: a ConfigParser object
        ConfigParser object for this project

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # check config
    #
    if config is None:
        logger.critical("Error: Need to supply a config")
        raise ValueError("Config is None")
    #
    # Set channel parameters
    #
    if spwtype == 'cont':
        restfreqs = [None for spw in my_spws.split(',')]
        start = None
        width = None
        nchan = None
        outframe = cp['contoutframe']
        veltype = None
        interpolation = None
    elif spwtype == 'line':
        spw_inds = [config.get("Spectral Windows","Line").split(',').index(spw) for spw in my_spws.split(',')]
        restfreqs = [config.get("Clean","restfreqs").split(',')[spw_ind] for spw_ind in spw_inds]
        start = cp['velstart']
        width = cp['chanwidth']
        nchan = cp['nchan']
        outframe = cp['lineoutframe']
        veltype = cp['veltype']
        interpolation = cp['interpolation']
    else:
        logger.critical("Error: spwtype {0} not supported".format(spwtype))
        raise ValueError("Invalid spwtype")
    #
    # Loop over spws
    #
    for spw,restfreq in zip(my_spws.split(','),restfreqs):
        #
        # Lightly clean spw
        #
        if regrid:
            myvis = vis.replace('.ms','')+'.spw{0}.cvel.ms'.format(spw)
            myspw = '0'
        else:
            myvis = vis
            myspw = spw
        imagename='{0}.spw{1}.channel.lightclean'.format(field,spw)
        if uvtaper:
            imagename = imagename + '.uvtaper'
            mask = '{0}.spw{1}.mfs.clean.uvtaper.mask'.format(field,spw)
        else:
            mask = '{0}.spw{1}.mfs.clean.mask'.format(field,spw)
        if not os.path.isdir(mask):
            logger.critical("Error: {0} does not exist".format(mask))
            raise ValueError("{0} does not exist".format(mask))
        logger.info("Lightly cleaning spw {0} (restfreq: {1})...".format(spw,restfreq))
        logger.info("Using mask: {0}".format(mask))
        casa.tclean(vis=myvis, imagename=imagename, field=field, spw=myspw, 
                    specmode='cube', threshold='0mJy', niter=cp['lightniter']*cp['nchan'], 
                    mask=mask, 
                    deconvolver='multiscale', scales=cp['scales'], 
                    gain=cp['gain'], cyclefactor=cp['cyclefactor'], 
                    imsize=cp['imsize'], pblimit=-1.0, cell=cp['cell'], 
                    weighting=cp['weighting'], robust=cp['robust'], 
                    restfreq=restfreq, start=start, width=width, nchan=nchan, 
                    outframe=outframe, veltype=veltype, interpolation=interpolation, 
                    uvtaper=cp['outertaper'], pbcor=False)
        #
        # Get RMS of residuals
        #
        dat = casa.imstat(imagename='{0}.residual'.format(imagename), 
                          axes=[0, 1], mask="'{0}.mask' == 0".format(imagename))
        logger.info("Max un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.max(dat['rms'])))
        logger.info("Mean un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.mean(dat['rms'])))
        logger.info("Median un-masked RMS: {0:.2f} mJy/beam".format(1000.*np.median(dat['rms'])))
        logger.info("Max un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.max(dat['medabsdevmed'])))
        logger.info("Mean un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.mean(dat['medabsdevmed'])))
        logger.info("Median un-masked MAD*1.4826: {0:.2f} mJy/beam".format(1000.*1.4826*np.median(dat['medabsdevmed'])))
        logger.info("Using median MADS*1.4826 times {0} (user-defined) as threshold".format(cp['nrms']))
        threshold = '{0:.2f}mJy'.format(cp["nrms"]*1000.*1.4826*np.median(dat['medabsdevmed']))
        #
        # Deep clean to threshold
        #
        imagename='{0}.spw{1}.channel.clean'.format(field, spw)
        if uvtaper:
            imagename = imagename + '.uvtaper'
            mask = '{0}.spw{1}.mfs.clean.uvtaper.mask'.format(field, spw)
        else:
            mask = '{0}.spw{1}.mfs.clean.mask'.format(field, spw)
        if not os.path.isdir(mask):
            logger.critical("Error: {0} does not exist".format(mask))
            raise ValueError("{0} does not exist".format(mask))
        logger.info("Cleaning spw {0} (restfreq: {1}) to threshold: {2}...".format(spw, restfreq, threshold))
        logger.info("Using mask: {0}".format(mask))
        casa.tclean(vis=myvis, imagename=imagename, field=field, spw=myspw, 
                    specmode='cube', threshold=threshold, niter=cp['maxniter']*cp['nchan'], 
                    mask=mask, 
                    deconvolver='multiscale', scales=cp['scales'], 
                    gain=cp['gain'], cyclefactor=cp['cyclefactor'], 
                    imsize=cp['imsize'], pblimit=-1.0, cell=cp['cell'], 
                    weighting=cp['weighting'], robust=cp['robust'], 
                    restfreq=restfreq, start=start, width=width, nchan=nchan, 
                    outframe=outframe, veltype=veltype, interpolation=interpolation, 
                    uvtaper=cp['outertaper'], pbcor=False, interactive=interactive)
        logger.info("Done.")
        #
        # Primary beam correction
        #
        pbimagename='{0}.spw{1}.channel.pb'.format(field, spw)
        if uvtaper:
            pbimagename = pbimagename + '.uvtaper'
        logger.info("Performing primary beam correction...")
        casa.impbcor(imagename='{0}.image'.format(imagename), 
                     pbimage='{0}.pb'.format(pbimagename), 
                     outfile='{0}.pbcor.image'.format(imagename))
        logger.info("Done.")
        #
        # Export to fits
        #
        logger.info("Exporting fits file...")
        casa.exportfits(imagename='{0}.image'.format(imagename), 
                        fitsimage='{0}.image.fits'.format(imagename), 
                        velocity=True, overwrite=True, history=False)
        casa.exportfits(imagename='{0}.mask'.format(imagename), 
                        fitsimage='{0}.mask.fits'.format(imagename), 
                        overwrite=True, history=False)
        casa.exportfits(imagename='{0}.residual'.format(imagename), 
                        fitsimage='{0}.residual.fits'.format(imagename), 
                        velocity=True, overwrite=True, history=False)
        casa.exportfits(imagename='{0}.pbcor.image'.format(imagename), 
                        fitsimage='{0}.pbcor.image.fits'.format(imagename), 
                        velocity=True, overwrite=True, history=False)
        logger.info("Done.")

def contplot(field, spws='', uvtaper=False, cp={}):
    """
    Generate PDF of MFS diagnostic plots

    Inputs:
      field :: string
        The field we're plotting
      spws :: string
        comma-separated string of spectral windows to plot
      uvtaper :: boolean
        if True, use uv-tapered images
      cp :: dictionary
        dictionary of clean parameters

    Returns:
      Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    logger.info("Generating continuum images...")
    #
    # Get center pixels
    #
    center_x = int(cp['imsize'][0]/2)
    center_y = int(cp['imsize'][1]/2)
    #
    # Loop over all plot filenames
    #
    goodplots = []
    for spw in spws.split(','):
        if spw != 'cont':
            spw = 'spw{0}'.format(spw)
        # check that this spectral window exists
        fname = '{0}.{1}.mfs.clean.image.fits'.format(field,spw)
        if uvtaper:
            fname = '{0}.{1}.mfs.clean.uvtaper.image.fits'.format(field,spw)
        if not os.path.exists(fname):
            continue
        if uvtaper:
            fitsfiles = ['{0}.{1}.mfs.dirty.uvtaper.image.fits'.format(field,spw),
                         '{0}.{1}.mfs.clean.uvtaper.image.fits'.format(field,spw),
                         '{0}.{1}.mfs.clean.uvtaper.residual.fits'.format(field,spw),
                         '{0}.{1}.mfs.clean.uvtaper.pbcor.image.fits'.format(field,spw)]
            maskfiles = ['{0}.{1}.mfs.clean.uvtaper.mask.fits'.format(field,spw),
                         '{0}.{1}.mfs.clean.uvtaper.mask.fits'.format(field,spw),
                         '{0}.{1}.mfs.clean.uvtaper.mask.fits'.format(field,spw),
                         '{0}.{1}.mfs.clean.uvtaper.mask.fits'.format(field,spw)]
            titles = ['{0} - spw: {1} - UVTap/Dirty'.format(field,spw),
                      '{0} - spw: {1} - UVTap/Clean'.format(field,spw),
                      '{0} - spw: {1} - UVTap/Residual'.format(field,spw),
                      '{0} - spw: {1} - UVTap/PBCorr'.format(field,spw)]
            labels = ['Flux Density (Jy/beam)',
                      'Flux Density (Jy/beam)',
                      'Flux Density (Jy/beam)',
                      'Flux Density (Jy/beam)']
            vlims = [(None,None),
                     (None,None),
                     (None,None),
                     (None,None)]
        else:
            fitsfiles = ['{0}.{1}.mfs.dirty.image.fits'.format(field,spw),
                         '{0}.{1}.mfs.clean.image.fits'.format(field,spw),
                         '{0}.{1}.mfs.clean.residual.fits'.format(field,spw),
                         '{0}.{1}.mfs.clean.pbcor.image.fits'.format(field,spw)]
            maskfiles = ['{0}.{1}.mfs.clean.mask.fits'.format(field,spw),
                         '{0}.{1}.mfs.clean.mask.fits'.format(field,spw),
                         '{0}.{1}.mfs.clean.mask.fits'.format(field,spw),
                         '{0}.{1}.mfs.clean.mask.fits'.format(field,spw)]
            titles = ['{0} - spw: {1} - Dirty'.format(field,spw),
                      '{0} - spw: {1} - Clean'.format(field,spw),
                      '{0} - spw: {1} - Residual'.format(field,spw),
                      '{0} - spw: {1} - PBCorr'.format(field,spw)]
            labels = ['Flux Density (Jy/beam)',
                      'Flux Density (Jy/beam)',
                      'Flux Density (Jy/beam)',
                      'Flux Density (Jy/beam)']
            vlims = [(None,None),
                     (None,None),
                     (None,None),
                     (None,None)]
        for fitsfile,maskfile,title,label,vlim in \
            zip(fitsfiles,maskfiles,titles,labels,vlims):
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
            cax = ax.imshow(hdu.data[0,0],
                            origin='lower',interpolation='none',
                            cmap='binary',vmin=vlim[0],vmax=vlim[1])
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
                cell = float(cp['cell'].replace('arcsec',''))
                beam_maj = hdu.header['BMAJ']*3600./cell
                beam_min = hdu.header['BMIN']*3600./cell
                beam_pa = hdu.header['BPA']
                ellipse = Ellipse((center_x-int(3.*center_x/4),
                                   center_y-int(3.*center_y/4)),
                                  beam_min,beam_maj,angle=beam_pa,
                                  fill=True,zorder=10,hatch='///',
                                  edgecolor='black',facecolor='white')
                ax.add_patch(ellipse)
            #
            # Plot colorbar
            #
            cbar = fig.colorbar(cax,fraction=0.046,pad=0.04)
            cbar.set_label(label)
            #
            # Re-scale to fit, then save
            #
            fname = fitsfile.replace('.fits','.pdf')
            fig.savefig(fname,bbox_inches='tight')
            plt.close(fig)
            plt.ion()
            goodplots.append(fname)
    #
    # Generate PDF of plots
    #
    # need to fix filenames so LaTeX doesn't complain
    outplots = ['{'+fn.replace('.pdf','')+'}.pdf' for fn in goodplots]
    logger.info("Generating PDF...")
    fname = '{0}.contplots.tex'.format(field)
    if uvtaper:
        fname = '{0}.uvtaper.contplots.tex'.format(field)
    with open(fname,'w') as f:
        f.write(r"\documentclass{article}"+"\n")
        f.write(r"\usepackage{graphicx}"+"\n")
        f.write(r"\usepackage[margin=0.1cm]{geometry}"+"\n")
        f.write(r"\begin{document}"+"\n")
        for i in range(0,len(outplots),4):
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
    logger.info("Done.")

def lineplot(field, line_spws='', uvtaper=False, cp={}):
    """
    Generate PDF of channel cube diagnostic plots

    Inputs:
      field :: string
        The field we're plotting
      line_spws :: string
        comma-separated string of spectral windows to plot
      uvtaper :: boolean
        if True, use uv-tapered images
      cp :: dictionary
        dictionary of clean parameters

    Returns:
      Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    logger.info("Generating line images...")
    #
    # Setup parameters for zooming
    #
    center_x = int(cp['imsize'][0]/2)
    center_y = int(cp['imsize'][1]/2)
    #
    # Get center pixels
    #
    goodplots = []
    for spw in line_spws.split(','):
        # check that this spectral window exists
        fname = '{0}.spw{1}.channel.clean.image.fits'.format(field,spw)
        if uvtaper:
            fname = '{0}.spw{1}.channel.clean.uvtaper.image.fits'.format(field,spw)
        if not os.path.exists(fname):
            continue
        #
        # Loop over all plot filenames
        #
        if uvtaper:
            fitsfiles = ['{0}.spw{1}.channel.dirty.uvtaper.image.fits'.format(field,spw),
                         '{0}.spw{1}.channel.clean.uvtaper.image.fits'.format(field,spw),
                         '{0}.spw{1}.channel.clean.uvtaper.residual.fits'.format(field,spw),
                         '{0}.spw{1}.channel.clean.uvtaper.pbcor.image.fits'.format(field,spw)]
            maskfile = '{0}.spw{1}.channel.clean.uvtaper.mask.fits'.format(field,spw)
            titles = ['{0} - spw: {1} - Dirty/UVTap'.format(field,spw),
                      '{0} - spw: {1} - Clean/UVTap'.format(field,spw),
                      '{0} - spw: {1} - Residual/UVTap'.format(field,spw),
                      '{0} - spw: {1} - UVTap/PBCorr'.format(field,spw)]
            labels = ['Flux Density (Jy/beam)',
                      'Flux Density (Jy/beam)',
                      'Flux Density (Jy/beam)',
                      'Flux Density (Jy/beam)']
            vlims = [(None,None),
                     (None,None),
                     (None,None),
                     (None,None)]
        else:
            fitsfiles = ['{0}.spw{1}.channel.dirty.image.fits'.format(field,spw),
                         '{0}.spw{1}.channel.clean.image.fits'.format(field,spw),
                         '{0}.spw{1}.channel.clean.residual.fits'.format(field,spw),
                         '{0}.spw{1}.channel.clean.pbcor.image.fits'.format(field,spw)]
            maskfile = '{0}.spw{1}.channel.clean.mask.fits'.format(field,spw)
            titles = ['{0} - spw: {1} - Dirty'.format(field,spw),
                      '{0} - spw: {1} - Clean'.format(field,spw),
                      '{0} - spw: {1} - Residual'.format(field,spw),
                      '{0} - spw: {1} - PBCorr'.format(field,spw)]
            labels = ['Flux Density (Jy/beam)',
                      'Flux Density (Jy/beam)',
                      'Flux Density (Jy/beam)',
                      'Flux Density (Jy/beam)']
            vlims = [(None,None),
                     (None,None),
                     (None,None),
                     (None,None)]
        for fitsfile,title,label,vlim in zip(fitsfiles,titles,labels,vlims):
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
                            origin='lower',interpolation='none',
                            cmap='binary',vmin=vlim[0],vmax=vlim[1])
            # ax.grid(True,color='black',ls='solid')
            ax.coords[0].set_major_formatter('hh:mm:ss')
            ax.set_xlabel('RA (J2000)')
            ax.set_ylabel('Declination (J2000)')
            #
            # Plot clean mask as contour
            #
            maskhdu = fits.open(maskfile)[0]
            contour = ax.contour(maskhdu.data[0,center_chan],levels=[0.5],
                                 origin='lower',colors='r',linewidths=2)
            #
            # Plot beam, if it is defined
            #
            if 'BMAJ' in hdu.header.keys():
                cell = float(cp['cell'].replace('arcsec',''))
                beam_maj = hdu.header['BMAJ']*3600./cell
                beam_min = hdu.header['BMIN']*3600./cell
                beam_pa = hdu.header['BPA']
                ellipse = Ellipse((center_x-int(3.*center_x/4),
                                center_y-int(3.*center_y/4)),
                                beam_min,beam_maj,angle=beam_pa,
                                fill=True,zorder=10,hatch='///',
                                edgecolor='black',facecolor='white')
                ax.add_patch(ellipse)
            elif len(hdulist) > 1:
                hdu = hdulist[1]
                cell = float(cp['cell'].replace('arcsec',''))
                beam_maj = hdu.data['BMAJ'][center_chan]/cell
                beam_min = hdu.data['BMIN'][center_chan]/cell
                beam_pa = hdu.data['BPA'][center_chan]
                ellipse = Ellipse((center_x-int(3.*center_x/4),
                                  center_y-int(3.*center_y/4)),
                                  beam_min,beam_maj,angle=beam_pa,
                                  fill=True,zorder=10,hatch='///',
                                  edgecolor='black',facecolor='white')
                ax.add_patch(ellipse)
            #
            # Plot colorbar
            #
            cbar = fig.colorbar(img,fraction=0.046,pad=0.04)
            cbar.set_label(label)
            #
            # Re-scale to fit, then save
            #
            fig.savefig(fitsfile.replace('.fits','.pdf'),
                        bbox_inches='tight')
            plt.close(fig)
            plt.ion()
            goodplots.append(fitsfile.replace('.fits','.pdf'))
        #
        # Generate spectrum
        #
        if uvtaper:
            fitsfile = '{0}.spw{1}.channel.clean.uvtaper.image.fits'.format(field,spw)
        else:
            fitsfile = '{0}.spw{1}.channel.clean.image.fits'.format(field,spw)
        hdu = fits.open(fitsfile)[0]
        spec = hdu.data[0,:,center_x,center_y]
        isnan = spec == 0.
        spec[isnan] = np.nan
        velo = (np.arange(len(spec))*hdu.header['CDELT3'] + hdu.header['CRVAL3'])/1000.
        #
        # Generate figure
        #
        plt.ioff()
        fig = plt.figure()
        ax = plt.subplot()
        ax.plot(velo,spec,'k-')
        ax.set_xlabel('Velocity (km/s)')
        ax.set_ylabel('Flux Density (Jy/beam)')
        ax.set_xlim(np.nanmin(velo),np.nanmax(velo))
        ybuff = 0.1*(np.nanmax(spec)-np.nanmin(spec))
        ax.set_ylim(np.nanmin(spec)-ybuff,np.nanmax(spec)+ybuff)
        if uvtaper:
            ax.set_title('{0} - spw {1} - UVTap/Center'.format(field,spw))
        else:
            ax.set_title('{0} - spw {1} - Center'.format(field,spw))
        ax.grid(False)
        fig.tight_layout()
        fig.savefig(fitsfile.replace('.fits','.spec.pdf'),
                    bbox_inches='tight')
        plt.close(fig)
        plt.ion()
        goodplots.append(fitsfile.replace('.fits','.spec.pdf'))
    #
    # Generate PDF of plots
    #
    # need to fix filenames so LaTeX doesn't complain
    goodplots = ['{'+fn.replace('.pdf','')+'}.pdf' for fn in goodplots]
    logger.info("Generating PDF...")
    fname = '{0}.lineplots.tex'.format(field)
    if uvtaper:
        fname = '{0}.uvtaper.lineplots.tex'.format(field)
    with open(fname,'w') as f:
        f.write(r"\documentclass{article}"+"\n")
        f.write(r"\usepackage{graphicx}"+"\n")
        f.write(r"\usepackage[margin=0.1cm]{geometry}"+"\n")
        f.write(r"\begin{document}"+"\n")
        for i in range(0,len(goodplots),5):
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
    logger.info("Done.")

def main(field, vis='', spws='', config_file='',
         uvtaper=False, regrid=False, auto='', skip_check=False,
         interactive=False):
    """
    Generate and clean images

    Inputs:
      field :: string
        The field name to clean
      vis :: string
        The measurement set containing all data for field
      spws :: string
        comma-separated list of line spws to clean
        if empty, clean all line spws
      config_file :: string
        filename of the configuration file for this project
      uvtaper :: boolean
        if True, apply UV tapering
      regrid :: boolean
        if True, use CVEL2 to re-grid the velocity axes
      skip_check :: boolean
        if True, skip checking for good spws
      interactive :: boolean
        if True, interactively clean
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
    if not os.path.exists(config_file):
        logger.critical('Configuration file not found')
        raise ValueError('Configuration file not found!')
    #
    # load configuration file
    #
    config = ConfigParser.ConfigParser()
    logger.info("Reading configuration file {0}".format(config_file))
    config.read(config_file)
    logger.info("Done.")
    #
    # initial setup, get all line and cont spws
    #
    my_cont_spws,all_line_spws,cp = setup(vis=vis,config=config,uvtaper=uvtaper)
    if spws == '':
        my_line_spws = all_line_spws
    else:
        my_line_spws = spws
    #
    # Determine which spws actually have data
    #
    if not skip_check:
        logger.info("Checking cont spws...")
        good_cont_spws = []
        for spw in my_cont_spws.split(','):
            foo = None
            foo = casa.visstat(vis=vis,spw=spw)
            if foo is not None:
                good_cont_spws.append(spw)
        my_cont_spws = ','.join(good_cont_spws)
        logger.info("Using cont spws: {0}".format(my_cont_spws))
        logger.info("Checking line spws...")
        good_line_spws = []
        for spw in my_line_spws.split(','):
            foo = None
            foo = casa.visstat(vis=vis,spw=spw)
            if foo is not None:
                good_line_spws.append(spw)
        my_line_spws = ','.join(good_line_spws)
        logger.info("Using line spws: {0}".format(my_line_spws))
    #
    # Re-grid
    #
    if regrid:
        regrid_velocity(vis=vis,spws=my_line_spws,config=config,cp=cp)
    #
    # Prompt the user with a menu for each option, or auto-do them
    #
    auto_items = auto.split(',')
    auto_ind = 0
    while True:
        if len(auto) == 0:
            print("0. Dirty image combined continuum spws (MFS; multi-term; multi-scale)")
            print("1. Autoclean combined continuum spws (MFS; multi-term; multi-scale)")
            print("2. Dirty image each continuum spw (MFS; multi-scale)")
            print("3. Autoclean each continuum spw (MFS; multi-scale)")
            print("4. Dirty image each continuum spw (channel; multi-scale)")
            print("5. Autoclean each continuum spw (channel; multi-scale)")
            print("6. Dirty image each line spw (MFS; multi-scale)")
            print("7. Autoclean each line spw (MFS; multi-scale)")
            print("8. Dirty image each line spw (channel; multi-scale)")
            print("9. Autoclean each line spw (channel; multi-scale)")
            print("10. Generate continuum and line diagnostic plots")
            print("q [quit]")
            answer = raw_input("> ")
        else:
            answer = auto_items[auto_ind]
            auto_ind += 1
        if answer == '0':
            mfs_dirty_cont(field=field, vis=vis, my_cont_spws=my_cont_spws, 
                           cp=cp, uvtaper=uvtaper)
        elif answer == '1':
            mfs_clean_cont(field=field, vis=vis, my_cont_spws=my_cont_spws, 
                           cp=cp, uvtaper=uvtaper, interactive=interactive)
        elif answer == '2':
            mfs_dirty_spws(field=field, vis=vis, my_spws=my_cont_spws, 
                           cp=cp, uvtaper=uvtaper)
        elif answer == '3':
            mfs_clean_spws(field=field, vis=vis, my_spws=my_cont_spws, 
                           cp=cp, uvtaper=uvtaper, interactive=interactive)
        elif answer == '4':
            channel_dirty_spws(field=field, vis=vis, my_spws=my_cont_spws, 
                               cp=cp, spwtype='cont', regrid=False, 
                               config=config, uvtaper=uvtaper)
        elif answer == '5':
            channel_clean_spws(field=field, vis=vis, my_spws=my_cont_spws, 
                               cp=cp, spwtype='cont', regrid=False, 
                               config=config, uvtaper=uvtaper, 
                               interactive=interactive)
        elif answer == '6':
            mfs_dirty_spws(field=field, vis=vis, my_spws=my_line_spws, 
                           cp=cp, uvtaper=uvtaper)
        elif answer == '7':
            mfs_clean_spws(field=field, vis=vis, my_spws=my_line_spws, 
                           cp=cp, uvtaper=uvtaper, interactive=interactive)
        elif answer == '8':
            channel_dirty_spws(field=field, vis=vis, my_spws=my_line_spws, 
                               cp=cp, spwtype='line', regrid=regrid, 
                               config=config, uvtaper=uvtaper)
        elif answer == '9':
            channel_clean_spws(field=field, vis=vis, my_spws=my_line_spws, 
                               cp=cp, spwtype='line', regrid=regrid, 
                               config=config, uvtaper=uvtaper, 
                               interactive=interactive)
        elif answer == '10':
            all_spws = ['cont']+natural_sort(my_cont_spws.split(',')+my_line_spws.split(','))
            contplot(field,spws=','.join(all_spws),uvtaper=uvtaper,cp=cp)
            lineplot(field,line_spws=all_line_spws,uvtaper=uvtaper,cp=cp)
        elif answer.lower() == 'q' or answer.lower() == 'quit':
            break
        else:
            print("Input not recognized.")
        if auto_ind >= len(auto_items):
            break
