"""
channel_images.py - Generate data cubes

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
import numpy as np
from astropy.io import fits

from .utils import makePB

import __main__ as casa


def grid_parameters(img, spws):
    """
    Determine the velocity grid parameters. Options:
    1. Define "start", "width", and "nchan" in the config file
    2. Define "start", "end", and "width" in the config file
    3. Define "start" and "end" in the config file. Use native width

    Inputs:
        img :: Imaging object
            The Imaging object
        spws :: list
            spectral windows to image

    Returns: start, width, nchans
        start :: string
            Starting velocity
        width :: string
            Channel width
        nchans :: list of integers
            Number of channels for each spectral window
    """
    start = None
    if img.cp["start"]:
        start = "{0}km/s".format(img.cp["start"])
    width = None
    if img.cp["width"]:
        width = "{0}km/s".format(img.cp["width"])
    nchans = [None for spw in spws]
    if img.cp["nchan"]:
        nchans = [int(img.cp["nchan"]) for spw in spws]
    #
    # If start and end parameters are set, determine nchan
    #
    if img.cp["start"] and img.cp["end"]:
        #
        # Determine number of channels based on requested end
        # and width
        #
        if img.cp["width"]:
            nchans = [
                int(
                    abs(float(img.cp["end"]) - float(img.cp["start"]))
                    / float(img.cp["width"])
                )
                + 1
                for spw in spws
            ]
        else:
            #
            # Otherwise, get native velocity width from MS
            #
            nchans = []
            casa.msmd.open(img.vis)
            for spw in spws:
                center_hz = casa.msmd.meanfreq(int(spw))
                width_hz = np.mean(casa.msmd.chanwidths(int(spw)))
                width_kms = width_hz / center_hz * 299792.458  # km/s
                nchans.append(
                    int(abs(float(img.cp["end"]) - float(img.cp["start"])) / width_kms)
                    + 1
                )
            casa.msmd.close()
    return start, width, nchans


def channel_dirty_spws(img, spws, spwtype):
    """
    Dirty image each supplied spectral window (channel cube)

    Inputs:
        img :: Imaging object
            The Imaging object
        spws :: string
            comma-separated string of spws to image
        spwtype :: string
            'line' or 'cont', to determine which clean params to use

    Returns: Nothing
    """
    # Set channel parameters
    if spwtype == "cont":
        restfreqs = [None for spw in spws]
        start = None
        width = img.cp["contwidth"]
        nchans = [None for spw in spws]
        outframe = img.cp["contoutframe"]
        veltype = None
        interpolation = None
    elif spwtype == "line":
        spw_inds = [
            img.config.get("Spectral Windows", "Line").split(",").index(spw)
            for spw in spws
        ]
        restfreqs = [
            img.config.get("Clean", "restfreqs").split(",")[spw_ind]
            for spw_ind in spw_inds
        ]
        outframe = img.cp["lineoutframe"]
        veltype = img.cp["veltype"]
        interpolation = img.cp["interpolation"]

        # Determine velocity-gridding parameter
        start, width, nchans = grid_parameters(img, spws)
    else:
        img.logger.critical("Error: spwtype {0} not supported".format(spwtype))
        raise ValueError("Invalid spwtype")

    # Loop over spws
    for spw, restfreq, nchan in zip(spws, restfreqs, nchans):
        imagename = "{0}.spw{1}.{2}.channel".format(img.imfield, spw, img.stokes)
        if img.uvtaper:
            imagename = imagename + ".uvtaper"
        imagename = os.path.join(img.outdir, imagename)
        img.logger.info(
            "Dirty imaging spw {0} (restfreq: {1})...".format(spw, restfreq)
        )
        casa.tclean(
            vis=img.vis,
            imagename=imagename,
            phasecenter=img.cp["phasecenter"],
            field=img.field,
            spw=spw,
            specmode="cube",
            gridder=img.cp["gridder"],
            wprojplanes=img.cp["wprojplanes"],
            threshold="0mJy",
            niter=0,
            deconvolver="multiscale",
            scales=img.cp["scales"],
            gain=img.cp["gain"],
            cyclefactor=img.cp["cyclefactor"],
            imsize=img.cp["imsize"],
            pblimit=-1.0,
            cell=img.cp["cell"],
            weighting=img.cp["weighting"],
            robust=img.cp["robust"],
            restfreq=restfreq,
            start=start,
            width=width,
            nchan=nchan,
            outframe=outframe,
            veltype=veltype,
            interpolation=interpolation,
            uvtaper=img.outertaper,
            uvrange=img.uvrange,
            stokes=img.stokes,
            pbcor=False,
            parallel=img.parallel,
            chanchunks=img.cp["chanchunks"],
        )
        img.logger.info("Done.")

        # Generate primary beam image
        img.logger.info(
            "Generating primary beam image of spw {0} (channel)...".format(spw)
        )
        makePB(
            vis=img.vis,
            field=img.field,
            spw=spw,
            uvrange=img.uvrange,
            stokes=img.stokes,
            imtemplate="{0}.image".format(imagename),
            outimage="{0}.pb.image".format(imagename),
            pblimit=img.cp["pblimit"],
        )

        # Primary beam correction
        pbimage = "{0}.pb.image".format(imagename)
        img.logger.info("Performing primary beam correction...")
        if not os.path.exists(pbimage):
            raise ValueError(
                "Could not find {0}. Did you dirty MFS first?".format(pbimage)
            )
        casa.impbcor(
            imagename="{0}.image".format(imagename),
            pbimage=pbimage,
            outfile="{0}.pbcor.image".format(imagename),
            overwrite=True,
        )
        img.logger.info("Done.")

        # Export to fits
        img.logger.info("Exporting fits file...")
        velocity = spwtype == "line"
        casa.exportfits(
            imagename="{0}.pb.image".format(imagename),
            fitsimage="{0}.pb.fits".format(imagename),
            velocity=velocity,
            overwrite=True,
            history=False,
        )
        casa.exportfits(
            imagename="{0}.image".format(imagename),
            fitsimage="{0}.dirty.image.fits".format(imagename),
            velocity=velocity,
            overwrite=True,
            history=False,
        )
        casa.exportfits(
            imagename="{0}.residual".format(imagename),
            fitsimage="{0}.dirty.residual.fits".format(imagename),
            velocity=velocity,
            overwrite=True,
            history=False,
        )
        casa.exportfits(
            imagename="{0}.pbcor.image".format(imagename),
            fitsimage="{0}.dirty.pbcor.image.fits".format(imagename),
            velocity=velocity,
            overwrite=True,
            history=False,
        )
        img.logger.info("Done.")


def channel_clean_spws(img, spws, spwtype):
    """
    Clean all supplied spws by channel using clean mask from MFS
    images.

    Inputs:
        img :: Imaging object
            The imaging object
        spws :: string
            comma-separated string of spws to image
        spwtype :: string
            'line' or 'cont', to determine which clean params to use

    Returns: Nothing
    """
    #
    # Set channel parameters
    #
    if spwtype == "cont":
        restfreqs = [None for spw in spws]
        start = None
        width = img.cp["contwidth"]
        nchans = [None for spw in spws]
        outframe = img.cp["contoutframe"]
        veltype = None
        interpolation = None
    elif spwtype == "line":
        spw_inds = [
            img.config.get("Spectral Windows", "Line").split(",").index(spw)
            for spw in spws
        ]
        restfreqs = [
            img.config.get("Clean", "restfreqs").split(",")[spw_ind]
            for spw_ind in spw_inds
        ]
        outframe = img.cp["lineoutframe"]
        veltype = img.cp["veltype"]
        interpolation = img.cp["interpolation"]

        # Determine velocity-gridding parameter
        start, width, nchans = grid_parameters(img, spws)
    else:
        img.logger.critical("Error: spwtype {0} not supported".format(spwtype))
        raise ValueError("Invalid spwtype")

    # Loop over spws
    for spw, restfreq, nchan in zip(spws, restfreqs, nchans):
        # Get niters
        if nchan is None:
            # get number of channels from dirty image
            imagename = "{0}.spw{1}.{2}.channel".format(img.imfield, spw, img.stokes)
            if img.uvtaper:
                imagename = imagename + ".uvtaper"
            imagename = os.path.join(img.outdir, imagename)
            imagename = imagename + ".dirty.image.fits"
            if not os.path.exists(imagename):
                raise ValueError(
                    "Must create dirty channel cube first: {0}".format(imagename)
                )
            dirty_hdr = fits.getheader(imagename)
            lightniter = img.cp["lightniter"] * dirty_hdr["NAXIS3"] * len(img.stokes)
            niter = img.cp["maxniter"] * dirty_hdr["NAXIS3"] * len(img.stokes)
        else:
            lightniter = img.cp["lightniter"] * nchan * len(img.stokes)
            niter = img.cp["maxniter"] * nchan * len(img.stokes)

        # If not interactive, Lightly clean spw
        imagename = "{0}.spw{1}.{2}.channel".format(img.imfield, spw, img.stokes)
        if img.uvtaper:
            imagename = imagename + ".uvtaper"
            mask = "{0}.spw{1}.{2}.mfs.uvtaper.mask".format(
                img.imfield, spw, img.stokes
            )
        else:
            mask = "{0}.spw{1}.{2}.mfs.mask".format(img.imfield, spw, img.stokes)
        imagename = os.path.join(img.outdir, imagename)
        mask = os.path.join(img.outdir, mask)
        if not os.path.isdir(mask):
            img.logger.critical("Error: {0} does not exist".format(mask))
            raise ValueError("{0} does not exist".format(mask))
        if not img.interactive:
            img.logger.info(
                "Lightly cleaning spw {0} (restfreq: {1})...".format(spw, restfreq)
            )
            img.logger.info("Using mask: {0}".format(mask))
            casa.tclean(
                vis=img.vis,
                imagename=imagename,
                phasecenter=img.cp["phasecenter"],
                field=img.field,
                spw=spw,
                gridder=img.cp["gridder"],
                wprojplanes=img.cp["wprojplanes"],
                specmode="cube",
                threshold="0mJy",
                niter=lightniter,
                mask=mask,
                deconvolver="multiscale",
                scales=img.cp["scales"],
                gain=img.cp["gain"],
                cyclefactor=img.cp["cyclefactor"],
                imsize=img.cp["imsize"],
                pblimit=-1.0,
                cell=img.cp["cell"],
                weighting=img.cp["weighting"],
                robust=img.cp["robust"],
                restfreq=restfreq,
                start=start,
                width=width,
                nchan=nchan,
                outframe=outframe,
                veltype=veltype,
                interpolation=interpolation,
                uvtaper=img.outertaper,
                uvrange=img.uvrange,
                stokes=img.stokes,
                pbcor=False,
                restart=True,
                calcres=False,
                calcpsf=False,
                parallel=img.parallel,
                chanchunks=img.cp["chanchunks"],
            )

            # This generates a channel mask, so next clean can't
            # have mfs mask
            mask = ""

            # Get RMS of residuals
            dat = casa.imstat(
                imagename="{0}.residual".format(imagename),
                axes=[0, 1],
                mask="'{0}.mask' == 0".format(imagename),
            )
            img.logger.info(
                "Max un-masked RMS: {0:.2f} mJy/beam".format(
                    1000.0 * np.max(dat["rms"])
                )
            )
            img.logger.info(
                "Max un-masked MAD*1.4826: {0:.2f} mJy/beam".format(
                    1000.0 * 1.4826 * np.max(dat["medabsdevmed"])
                )
            )
            img.logger.info(
                "Using max(MAD) x 1.4826 x {0} as threshold".format(img.cp["nrms"])
            )
            threshold = "{0:.2f}mJy".format(
                img.cp["nrms"] * 1000.0 * 1.4826 * np.max(dat["medabsdevmed"])
            )
        else:
            threshold = "0.0mJy"

        # Deep clean to threshold
        img.logger.info(
            "Cleaning spw {0} (restfreq: {1}) to threshold: {2}...".format(
                spw, restfreq, threshold
            )
        )
        if mask:
            img.logger.info("Using mask: {0}".format(mask))
        casa.tclean(
            vis=img.vis,
            imagename=imagename,
            phasecenter=img.cp["phasecenter"],
            field=img.field,
            spw=spw,
            specmode="cube",
            gridder=img.cp["gridder"],
            wprojplanes=img.cp["wprojplanes"],
            threshold=threshold,
            niter=niter,
            mask=mask,
            deconvolver="multiscale",
            scales=img.cp["scales"],
            gain=img.cp["gain"],
            cyclefactor=img.cp["cyclefactor"],
            imsize=img.cp["imsize"],
            pblimit=-1.0,
            cell=img.cp["cell"],
            weighting=img.cp["weighting"],
            robust=img.cp["robust"],
            restfreq=restfreq,
            start=start,
            width=width,
            nchan=nchan,
            outframe=outframe,
            veltype=veltype,
            interpolation=interpolation,
            uvtaper=img.outertaper,
            uvrange=img.uvrange,
            pbcor=False,
            stokes=img.stokes,
            interactive=img.interactive,
            restart=True,
            calcres=False,
            calcpsf=False,
            parallel=img.parallel,
            chanchunks=img.cp["chanchunks"],
        )
        img.logger.info("Done.")

        # Primary beam correction
        pbimage = "{0}.pb.image".format(imagename)
        img.logger.info("Performing primary beam correction...")
        if not os.path.exists(pbimage):
            raise ValueError(
                "Could not find {0}. Did you dirty MFS first?".format(pbimage)
            )
        casa.impbcor(
            imagename="{0}.image".format(imagename),
            pbimage=pbimage,
            outfile="{0}.pbcor.image".format(imagename),
            overwrite=True,
        )
        img.logger.info("Done.")

        # Export to fits
        img.logger.info("Exporting fits file...")
        velocity = spwtype == "line"
        casa.exportfits(
            imagename="{0}.image".format(imagename),
            fitsimage="{0}.clean.image.fits".format(imagename),
            velocity=velocity,
            overwrite=True,
            history=False,
        )
        casa.exportfits(
            imagename="{0}.mask".format(imagename),
            fitsimage="{0}.clean.mask.fits".format(imagename),
            overwrite=True,
            history=False,
        )
        casa.exportfits(
            imagename="{0}.residual".format(imagename),
            fitsimage="{0}.clean.residual.fits".format(imagename),
            velocity=velocity,
            overwrite=True,
            history=False,
        )
        casa.exportfits(
            imagename="{0}.pbcor.image".format(imagename),
            fitsimage="{0}.clean.pbcor.image.fits".format(imagename),
            velocity=velocity,
            overwrite=True,
            history=False,
        )
        img.logger.info("Done.")
