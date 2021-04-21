"""
mfs_images.py - Generate multi-frequency synthesis images

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

from .utils import makePB

import __main__ as casa


def mfs_dirty_cont(img):
    """
    Dirty image combined spws using multi-frequency synthesis

    Inputs:
        img :: Imaging object
            The Imaging object

    Returns: Nothing
    """
    imagename = "{0}.cont.{1}.mfs".format(img.field, img.stokes)
    if img.uvtaper:
        imagename = imagename + ".uvtaper"
    img.logger.info("Generating dirty continuum image (MFS)...")
    cleanspw = ",".join(
        [
            "{0}:{1}".format(spw, chan) if chan else spw
            for spw, chan in zip(img.cont_spws, img.cont_chans)
        ]
    )
    casa.tclean(
        vis=img.vis,
        imagename=os.path.join(img.outdir, imagename),
        phasecenter=img.cp["phasecenter"],
        field=img.field,
        spw=cleanspw,
        gridder=img.cp["gridder"],
        wprojplanes=img.cp["wprojplanes"],
        specmode="mfs",
        threshold="0mJy",
        niter=0,
        nterms=img.cp["nterms"],
        deconvolver="mtmfs",
        scales=img.cp["scales"],
        gain=img.cp["gain"],
        cyclefactor=img.cp["cyclefactor"],
        imsize=img.cp["imsize"],
        pblimit=-1.0,
        cell=img.cp["cell"],
        weighting=img.cp["weighting"],
        robust=img.cp["robust"],
        uvtaper=img.outertaper,
        uvrange=img.uvrange,
        stokes=img.stokes,
        pbcor=False,
        parallel=img.parallel,
    )
    img.logger.info("Done.")

    # Primary beam correction
    img.logger.info("Performing primary beam correction...")
    # Due to widebandpbcor limitiations, need to go into outdir
    os.chdir(img.outdir)
    spwlist = [
        int(spw) for chan in img.cp["contpbchan"].split(",") for spw in img.cont_spws
    ]
    chanlist = [
        int(chan) for chan in img.cp["contpbchan"].split(",") for spw in img.cont_spws
    ]
    weightlist = [1.0 for _ in spwlist]
    casa.widebandpbcor(
        vis=os.path.join("..", img.vis),
        imagename=imagename,
        nterms=img.cp["nterms"],
        pbmin=img.cp["pblimit"],
        threshold="0.1mJy",
        spwlist=spwlist,
        weightlist=weightlist,
        chanlist=chanlist,
    )
    img.logger.info("Done.")
    os.chdir("..")

    # Export to fits
    img.logger.info("Exporting fits file...")
    casa.exportfits(
        imagename="{0}/{1}.pbcor.workdirectory/{1}.pb.tt0".format(
            img.outdir, imagename
        ),
        fitsimage="{0}/{1}.pb.fits".format(img.outdir, imagename),
        overwrite=True,
        history=False,
    )
    imagename = os.path.join(img.outdir, imagename)
    casa.exportfits(
        imagename="{0}.image.tt0".format(imagename),
        fitsimage="{0}.dirty.image.fits".format(imagename),
        overwrite=True,
        history=False,
    )
    casa.exportfits(
        imagename="{0}.residual.tt0".format(imagename),
        fitsimage="{0}.dirty.residual.fits".format(imagename),
        overwrite=True,
        history=False,
    )
    casa.exportfits(
        imagename="{0}.pbcor.image.tt0".format(imagename),
        fitsimage="{0}.dirty.pbcor.image.fits".format(imagename),
        overwrite=True,
        history=False,
    )
    casa.exportfits(
        imagename="{0}.pbcor.image.alpha".format(imagename),
        fitsimage="{0}.dirty.pbcor.image.alpha.fits".format(imagename),
        overwrite=True,
        history=False,
    )
    casa.exportfits(
        imagename="{0}.pbcor.image.alpha.error".format(imagename),
        fitsimage="{0}.dirty.pbcor.image.alpha.error.fits".format(imagename),
        overwrite=True,
        history=False,
    )
    img.logger.info("Done.")


def mfs_clean_cont(img):
    """
    Clean combined spws using multi-frequency synthesis

    Inputs:
        img :: Imaging object
            The Imaging object

    Returns: Nothing
    """
    # If not cleaning interactively, lightly clean to get RMS
    # threshold
    imagename = "{0}.cont.{1}.mfs".format(img.field, img.stokes)
    if img.uvtaper:
        imagename = imagename + ".uvtaper"
    cleanspw = ",".join(
        [
            "{0}:{1}".format(spw, chan) if chan else spw
            for spw, chan in zip(img.cont_spws, img.cont_chans)
        ]
    )
    if not img.interactive:
        savemodel = "none"
        if img.savemodel == "light":
            savemodel = "modelcolumn"
        img.logger.info("Lightly cleaning continuum image (MFS)...")
        casa.tclean(
            vis=img.vis,
            imagename=os.path.join(img.outdir, imagename),
            phasecenter=img.cp["phasecenter"],
            field=img.field,
            spw=cleanspw,
            gridder=img.cp["gridder"],
            wprojplanes=img.cp["wprojplanes"],
            specmode="mfs",
            threshold="0mJy",
            niter=img.cp["lightniter"] * len(img.stokes),
            usemask="auto-multithresh",
            pbmask=img.cp["contpbmask"],
            sidelobethreshold=img.cp["contsidelobethreshold"],
            noisethreshold=img.cp["contnoisethreshold"],
            lownoisethreshold=img.cp["contlownoisethreshold"],
            negativethreshold=img.cp["contnegativethreshold"],
            smoothfactor=img.cp["contsmoothfactor"],
            minbeamfrac=img.cp["contminbeamfrac"],
            cutthreshold=img.cp["contcutthreshold"],
            growiterations=img.cp["contgrowiterations"],
            nterms=img.cp["nterms"],
            deconvolver="mtmfs",
            scales=img.cp["scales"],
            gain=img.cp["gain"],
            cyclefactor=img.cp["cyclefactor"],
            imsize=img.cp["imsize"],
            pblimit=-1.0,
            cell=img.cp["cell"],
            weighting=img.cp["weighting"],
            robust=img.cp["robust"],
            uvtaper=img.outertaper,
            uvrange=img.uvrange,
            stokes=img.stokes,
            pbcor=False,
            restart=True,
            calcres=False,
            calcpsf=False,
            savemodel=savemodel,
            parallel=img.parallel,
        )
        img.logger.info("Done.")

        # Get RMS of residuals outside of clean mask
        dat = casa.imstat(
            imagename="{0}.residual.tt0".format(os.path.join(img.outdir, imagename)),
            axes=[0, 1],
            mask="'{0}.mask' == 0".format(os.path.join(img.outdir, imagename)),
        )
        img.logger.info(
            "Max un-masked RMS: {0:.2f} mJy/beam".format(1000.0 * np.max(dat["rms"]))
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
        # No threshold for interactive clean
        threshold = "0.0mJy"

    # Clean to threshold
    savemodel = "none"
    if img.savemodel == "clean":
        savemodel = "modelcolumn"
    img.logger.info(
        "Cleaning continuum image (MFS) to threshold: {0}...".format(threshold)
    )
    casa.tclean(
        vis=img.vis,
        imagename=os.path.join(img.outdir, imagename),
        phasecenter=img.cp["phasecenter"],
        field=img.field,
        spw=cleanspw,
        gridder=img.cp["gridder"],
        wprojplanes=img.cp["wprojplanes"],
        specmode="mfs",
        threshold=threshold,
        niter=img.cp["maxniter"] * len(img.stokes),
        usemask="auto-multithresh",
        pbmask=img.cp["contpbmask"],
        sidelobethreshold=img.cp["contsidelobethreshold"],
        noisethreshold=img.cp["contnoisethreshold"],
        lownoisethreshold=img.cp["contlownoisethreshold"],
        negativethreshold=img.cp["contnegativethreshold"],
        smoothfactor=img.cp["contsmoothfactor"],
        minbeamfrac=img.cp["contminbeamfrac"],
        cutthreshold=img.cp["contcutthreshold"],
        growiterations=img.cp["contgrowiterations"],
        nterms=img.cp["nterms"],
        deconvolver="mtmfs",
        scales=img.cp["scales"],
        gain=img.cp["gain"],
        cyclefactor=img.cp["cyclefactor"],
        imsize=img.cp["imsize"],
        pblimit=-1.0,
        cell=img.cp["cell"],
        weighting=img.cp["weighting"],
        robust=img.cp["robust"],
        uvtaper=img.outertaper,
        uvrange=img.uvrange,
        pbcor=False,
        stokes=img.stokes,
        interactive=img.interactive,
        restart=True,
        calcres=False,
        calcpsf=False,
        savemodel=savemodel,
        parallel=img.parallel,
    )
    img.logger.info("Done.")

    # Primary beam correction using PB of center channel
    img.logger.info("Performing primary beam correction...")
    # Due to widebandpbcor limitiations, need to go into outdir
    os.chdir(img.outdir)
    spwlist = [
        int(spw) for chan in img.cp["contpbchan"].split(",") for spw in img.cont_spws
    ]
    chanlist = [
        int(chan) for chan in img.cp["contpbchan"].split(",") for spw in img.cont_spws
    ]
    weightlist = [1.0 for _ in spwlist]
    casa.widebandpbcor(
        vis=os.path.join("..", img.vis),
        imagename=imagename,
        nterms=img.cp["nterms"],
        pbmin=img.cp["pblimit"],
        threshold="0.1mJy",
        spwlist=spwlist,
        weightlist=weightlist,
        chanlist=chanlist,
    )
    img.logger.info("Done.")
    os.chdir("..")

    # Export to fits
    img.logger.info("Exporting fits file...")
    imagename = os.path.join(img.outdir, imagename)
    casa.exportfits(
        imagename="{0}.image.tt0".format(imagename),
        fitsimage="{0}.clean.image.fits".format(imagename),
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
        imagename="{0}.residual.tt0".format(imagename),
        fitsimage="{0}.clean.residual.fits".format(imagename),
        overwrite=True,
        history=False,
    )
    casa.exportfits(
        imagename="{0}.pbcor.image.tt0".format(imagename),
        fitsimage="{0}.clean.pbcor.image.fits".format(imagename),
        overwrite=True,
        history=False,
    )
    casa.exportfits(
        imagename="{0}.pbcor.image.alpha".format(imagename),
        fitsimage="{0}.clean.pbcor.image.alpha.fits".format(imagename),
        overwrite=True,
        history=False,
    )
    casa.exportfits(
        imagename="{0}.pbcor.image.alpha.error".format(imagename),
        fitsimage="{0}.clean.pbcor.image.alpha.error.fits".format(imagename),
        overwrite=True,
        history=False,
    )
    img.logger.info("Done.")


def mfs_dirty_spws(img, spws, chans):
    """
    Dirty image each supplied spw (MFS)

    Inputs:
        img :: Imaging object
            The Imaging object
        spws, chans :: lists
            spectral windows and channels to image

    Returns: Nothing
    """
    for spw, chan in zip(spws, chans):
        imagename = "{0}.spw{1}.{2}.mfs".format(img.field, spw, img.stokes)
        if img.uvtaper:
            imagename = imagename + ".uvtaper"
        imagename = os.path.join(img.outdir, imagename)
        img.logger.info("Generating dirty image of spw {0} (MFS)...".format(spw))
        cleanspw = spw
        if chan:
            cleanspw = "{0}:{1}".format(spw, chan)
        casa.tclean(
            vis=img.vis,
            imagename=imagename,
            phasecenter=img.cp["phasecenter"],
            field=img.field,
            spw=cleanspw,
            specmode="mfs",
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
            uvtaper=img.outertaper,
            uvrange=img.uvrange,
            stokes=img.stokes,
            pbcor=False,
            parallel=img.parallel,
        )
        img.logger.info("Done.")

        # Generate primary beam image
        img.logger.info("Generating primary beam image of spw {0} (MFS)...".format(spw))
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
        # drop degenerate axes (channel and stokes if stokes=I)
        casa.imsubimage(
            imagename="{0}.pb.image".format(imagename),
            outfile="{0}.pb.image.sub".format(imagename),
            dropdeg=True,
            overwrite=True,
        )
        img.logger.info("Done.")

        # Primary beam correction
        img.logger.info("Performing primary beam correction...")
        casa.impbcor(
            imagename="{0}.image".format(imagename),
            pbimage="{0}.pb.image".format(imagename),
            outfile="{0}.pbcor.image".format(imagename),
            overwrite=True,
        )
        img.logger.info("Done.")

        # Export to fits
        img.logger.info("Exporting fits file...")
        casa.exportfits(
            imagename="{0}.pb.image.sub".format(imagename),
            fitsimage="{0}.pb.fits".format(imagename),
            overwrite=True,
            history=False,
        )
        casa.exportfits(
            imagename="{0}.image".format(imagename),
            fitsimage="{0}.dirty.image.fits".format(imagename),
            overwrite=True,
            history=False,
        )
        casa.exportfits(
            imagename="{0}.residual".format(imagename),
            fitsimage="{0}.dirty.residual.fits".format(imagename),
            overwrite=True,
            history=False,
        )
        casa.exportfits(
            imagename="{0}.pbcor.image".format(imagename),
            fitsimage="{0}.dirty.pbcor.image.fits".format(imagename),
            overwrite=True,
            history=False,
        )
        img.logger.info("Done.")


def mfs_clean_spws(img, spws, chans, spwtype):
    """
    Clean each supplied spw (MFS)

    Inputs:
        img :: Imaging object
            The Imaging object
        spws, chans :: lists
            spectral windows and channels to image
        spwtype :: string
            'line' or 'cont', to determine which clean params to use

    Returns: Nothing
    """
    for spw, chan in zip(spws, chans):
        # If not interactive, Lightly clean to get threshold
        imagename = "{0}.spw{1}.{2}.mfs".format(img.field, spw, img.stokes)
        if img.uvtaper:
            imagename = imagename + ".uvtaper"
        imagename = os.path.join(img.outdir, imagename)
        cleanspw = spw
        if chan:
            cleanspw = "{0}:{1}".format(spw, chan)
        if not img.interactive:
            # Save model if necessary
            savemodel = "none"
            if img.savemodel == "light":
                savemodel = "modelcolumn"
            img.logger.info("Lightly cleaning spw {0} (MFS)...".format(spw))
            casa.tclean(
                vis=img.vis,
                imagename=imagename,
                phasecenter=img.cp["phasecenter"],
                field=img.field,
                spw=cleanspw,
                specmode="mfs",
                gridder=img.cp["gridder"],
                wprojplanes=img.cp["wprojplanes"],
                threshold="0mJy",
                niter=img.cp["lightniter"] * len(img.stokes),
                usemask="auto-multithresh",
                pbmask=img.cp[spwtype + "pbmask"],
                sidelobethreshold=img.cp[spwtype + "sidelobethreshold"],
                noisethreshold=img.cp[spwtype + "noisethreshold"],
                lownoisethreshold=img.cp[spwtype + "lownoisethreshold"],
                negativethreshold=img.cp[spwtype + "negativethreshold"],
                smoothfactor=img.cp[spwtype + "smoothfactor"],
                minbeamfrac=img.cp[spwtype + "minbeamfrac"],
                cutthreshold=img.cp[spwtype + "cutthreshold"],
                growiterations=img.cp[spwtype + "growiterations"],
                deconvolver="multiscale",
                scales=img.cp["scales"],
                gain=img.cp["gain"],
                cyclefactor=img.cp["cyclefactor"],
                imsize=img.cp["imsize"],
                pblimit=-1.0,
                cell=img.cp["cell"],
                weighting=img.cp["weighting"],
                robust=img.cp["robust"],
                uvtaper=img.outertaper,
                uvrange=img.uvrange,
                stokes=img.stokes,
                savemodel=savemodel,
                pbcor=False,
                restart=True,
                calcres=False,
                calcpsf=False,
                parallel=img.parallel,
            )
            img.logger.info("Done.")

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

        # Clean to threshold
        # Save model if necessary
        savemodel = "none"
        if img.savemodel == "clean":
            savemodel = "modelcolumn"
        img.logger.info(
            "Cleaning spw {0} (MFS) to threshold: {1}...".format(spw, threshold)
        )
        casa.tclean(
            vis=img.vis,
            imagename=imagename,
            field=img.field,
            phasecenter=img.cp["phasecenter"],
            spw=cleanspw,
            gridder=img.cp["gridder"],
            wprojplanes=img.cp["wprojplanes"],
            specmode="mfs",
            threshold=threshold,
            niter=img.cp["maxniter"] * len(img.stokes),
            usemask="auto-multithresh",
            pbmask=img.cp[spwtype + "pbmask"],
            sidelobethreshold=img.cp[spwtype + "sidelobethreshold"],
            noisethreshold=img.cp[spwtype + "noisethreshold"],
            lownoisethreshold=img.cp[spwtype + "lownoisethreshold"],
            negativethreshold=img.cp[spwtype + "negativethreshold"],
            smoothfactor=img.cp[spwtype + "smoothfactor"],
            minbeamfrac=img.cp[spwtype + "minbeamfrac"],
            cutthreshold=img.cp[spwtype + "cutthreshold"],
            growiterations=img.cp[spwtype + "growiterations"],
            deconvolver="multiscale",
            scales=img.cp["scales"],
            gain=img.cp["gain"],
            cyclefactor=img.cp["cyclefactor"],
            imsize=img.cp["imsize"],
            pblimit=-1.0,
            cell=img.cp["cell"],
            weighting=img.cp["weighting"],
            robust=img.cp["robust"],
            uvtaper=img.outertaper,
            uvrange=img.uvrange,
            pbcor=False,
            stokes=img.stokes,
            savemodel=savemodel,
            interactive=img.interactive,
            restart=True,
            calcres=False,
            calcpsf=False,
            parallel=img.parallel,
        )
        img.logger.info("Done.")

        # Primary beam correction
        img.logger.info("Performing primary beam correction...")
        casa.impbcor(
            imagename="{0}.image".format(imagename),
            pbimage="{0}.pb.image".format(imagename),
            outfile="{0}.pbcor.image".format(imagename),
            overwrite=True,
        )
        img.logger.info("Done.")

        # Export to fits
        img.logger.info("Exporting fits file...")
        casa.exportfits(
            imagename="{0}.image".format(imagename),
            fitsimage="{0}.clean.image.fits".format(imagename),
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
            overwrite=True,
            history=False,
        )
        casa.exportfits(
            imagename="{0}.pbcor.image".format(imagename),
            fitsimage="{0}.clean.pbcor.image.fits".format(imagename),
            overwrite=True,
            history=False,
        )
        img.logger.info("Done.")
