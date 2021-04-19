"""
imaging_plots.py - WISP plotting utilities for imaging

Copyright(C) 2018-2020 by
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
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from astropy.io import fits
from astropy.wcs import WCS

from .utils import natural_sort


def contplot(img):
    """
    Generate PDF of MFS diagnostic plots

    Inputs:
        img :: Imaging object
            The imaging object

    Returns: Nothing
    """

    # Get center pixels
    center_x = int(img.cp["imsize"][0] / 2)
    center_y = int(img.cp["imsize"][1] / 2)

    # Loop over all plot filenames
    goodplots = []
    spws = ["cont"] + natural_sort(img.cont_spws + img.line_spws)
    for spw in spws:
        if spw != "cont":
            spw = "spw{0}".format(spw)

        # check that this spectral window exists
        fname = "{0}.{1}.{2}.mfs.clean.image.fits".format(img.field, spw, img.stokes)
        if img.uvtaper:
            fname = "{0}.{1}.{2}.mfs.uvtaper.clean.image.fits".format(
                img.field, spw, img.stokes
            )
        fname = os.path.join(img.outdir, fname)
        if not os.path.exists(fname):
            continue
        if img.uvtaper:
            fitsfiles = [
                "{0}.{1}.{2}.mfs.uvtaper.dirty.image.fits".format(
                    img.field, spw, img.stokes
                ),
                "{0}.{1}.{2}.mfs.uvtaper.clean.image.fits".format(
                    img.field, spw, img.stokes
                ),
                "{0}.{1}.{2}.mfs.uvtaper.clean.residual.fits".format(
                    img.field, spw, img.stokes
                ),
                "{0}.{1}.{2}.mfs.uvtaper.clean.pbcor.image.fits".format(
                    img.field, spw, img.stokes
                ),
            ]
            maskfiles = [
                "{0}.{1}.{2}.mfs.uvtaper.clean.mask.fits".format(
                    img.field, spw, img.stokes
                ),
                "{0}.{1}.{2}.mfs.uvtaper.clean.mask.fits".format(
                    img.field, spw, img.stokes
                ),
                "{0}.{1}.{2}.mfs.uvtaper.clean.mask.fits".format(
                    img.field, spw, img.stokes
                ),
                "{0}.{1}.{2}.mfs.uvtaper.clean.mask.fits".format(
                    img.field, spw, img.stokes
                ),
            ]
            titles = [
                "{0} - {1} - Taper/Dirty".format(img.field, spw),
                "{0} - {1} - Taper/Clean".format(img.field, spw),
                "{0} - {1} - Taper/Residual".format(img.field, spw),
                "{0} - {1} - Taper/PBCorr".format(img.field, spw),
            ]
            labels = [
                "Flux Density (Jy/beam)",
                "Flux Density (Jy/beam)",
                "Flux Density (Jy/beam)",
                "Flux Density (Jy/beam)",
            ]
            vlims = [
                (None, None),
                (None, None),
                (None, None),
                (None, None),
            ]
        else:
            fitsfiles = [
                "{0}.{1}.{2}.mfs.dirty.image.fits".format(img.field, spw, img.stokes),
                "{0}.{1}.{2}.mfs.clean.image.fits".format(img.field, spw, img.stokes),
                "{0}.{1}.{2}.mfs.clean.residual.fits".format(
                    img.field, spw, img.stokes
                ),
                "{0}.{1}.{2}.mfs.clean.pbcor.image.fits".format(
                    img.field, spw, img.stokes
                ),
            ]
            maskfiles = [
                "{0}.{1}.{2}.mfs.clean.mask.fits".format(img.field, spw, img.stokes),
                "{0}.{1}.{2}.mfs.clean.mask.fits".format(img.field, spw, img.stokes),
                "{0}.{1}.{2}.mfs.clean.mask.fits".format(img.field, spw, img.stokes),
                "{0}.{1}.{2}.mfs.clean.mask.fits".format(img.field, spw, img.stokes),
            ]
            titles = [
                "{0} - {1} - Dirty".format(img.field, spw),
                "{0} - {1} - Clean".format(img.field, spw),
                "{0} - {1} - Residual".format(img.field, spw),
                "{0} - {1} - PBCorr".format(img.field, spw),
            ]
            labels = [
                "Flux Density (Jy/beam)",
                "Flux Density (Jy/beam)",
                "Flux Density (Jy/beam)",
                "Flux Density (Jy/beam)",
            ]
            vlims = [
                (None, None),
                (None, None),
                (None, None),
                (None, None),
            ]

        # Generate figures for each Stokes parameter
        for ind, stokes in enumerate(img.stokes):

            # Loop over figures
            for fitsfile, maskfile, title, label, vlim in zip(
                fitsfiles, maskfiles, titles, labels, vlims
            ):
                # Open fits file, generate WCS
                hdulist = fits.open(os.path.join(img.outdir, fitsfile))
                hdu = hdulist[0]
                wcs = WCS(hdu.header)

                # Set-up figure
                plt.ioff()
                fig = plt.figure()
                ax = plt.subplot(projection=wcs.sub(["celestial"]))
                ax.set_title("{0} - {1}".format(title, stokes))
                cax = ax.imshow(
                    hdu.data[ind, 0],
                    origin="lower",
                    interpolation="none",
                    cmap="binary",
                    vmin=vlim[0],
                    vmax=vlim[1],
                )
                # ax.grid(True,color='black',ls='solid')
                if img.cp["frame"] == "J2000":
                    ax.coords[0].set_major_formatter("hh:mm:ss")
                    ax.set_xlabel("RA (J2000)")
                    ax.set_ylabel("Declination (J2000)")
                elif img.cp["frame"] == "GALACTIC":
                    ax.set_xlabel("Galactic Longitude")
                    ax.set_ylabel("Galactic Latitude")

                # Plot clean mask as contour
                if maskfile != "":
                    maskhdu = fits.open(os.path.join(img.outdir, maskfile))[0]
                    _ = ax.contour(
                        maskhdu.data[ind, 0],
                        levels=[0.5],
                        origin="lower",
                        colors="r",
                        linewidths=2,
                    )

                # Plot beam, if it is defined
                if "BMAJ" in hdu.header.keys():
                    cell = float(img.cp["cell"].replace("arcsec", ""))
                    beam_maj = hdu.header["BMAJ"] * 3600.0 / cell
                    beam_min = hdu.header["BMIN"] * 3600.0 / cell
                    beam_pa = hdu.header["BPA"]
                    ellipse = Ellipse(
                        (
                            center_x - int(3.0 * center_x / 4),
                            center_y - int(3.0 * center_y / 4),
                        ),
                        beam_min,
                        beam_maj,
                        angle=beam_pa,
                        fill=True,
                        zorder=10,
                        hatch="///",
                        edgecolor="black",
                        facecolor="white",
                    )
                    ax.add_patch(ellipse)
                elif len(hdulist) > 1:
                    hdu = hdulist[1]
                    cell = float(img.cp["cell"].replace("arcsec", ""))
                    beam_maj = hdu.data["BMAJ"][ind] / cell
                    beam_min = hdu.data["BMIN"][ind] / cell
                    beam_pa = hdu.data["BPA"][ind]
                    ellipse = Ellipse(
                        (
                            center_x - int(3.0 * center_x / 4),
                            center_y - int(3.0 * center_y / 4),
                        ),
                        beam_min,
                        beam_maj,
                        angle=beam_pa,
                        fill=True,
                        zorder=10,
                        hatch="///",
                        edgecolor="black",
                        facecolor="white",
                    )
                    ax.add_patch(ellipse)

                # Plot colorbar
                cbar = fig.colorbar(cax, fraction=0.046, pad=0.04)
                cbar.set_label(label)

                # Re-scale to fit, then save
                fname = fitsfile.replace(".fits", ".{0}.pdf".format(stokes))
                fig.savefig(os.path.join(img.outdir, fname), bbox_inches="tight")
                plt.close(fig)
                plt.ion()
                goodplots.append(fname)
    #
    # Generate PDF of plots
    #
    # need to fix filenames so LaTeX doesn't complain
    outplots = ["{" + fn.replace(".pdf", "") + "}.pdf" for fn in goodplots]
    img.logger.info("Generating PDF...")
    fname = "{0}.contplots.tex".format(img.field)
    if img.uvtaper:
        fname = "{0}.uvtaper.contplots.tex".format(img.field)
    with open(os.path.join(img.outdir, fname), "w") as f:
        f.write(r"\documentclass{article}" + "\n")
        f.write(r"\usepackage{graphicx}" + "\n")
        f.write(r"\usepackage[margin=0.1cm]{geometry}" + "\n")
        f.write(r"\begin{document}" + "\n")
        for i in range(0, len(outplots), 4):
            f.write(r"\begin{figure}" + "\n")
            f.write(r"\centering" + "\n")
            f.write(
                r"\includegraphics[width=0.45\textwidth]{" + outplots[i] + r"}" + "\n"
            )
            f.write(
                r"\includegraphics[width=0.45\textwidth]{"
                + outplots[i + 1]
                + r"} \\"
                + "\n"
            )
            f.write(
                r"\includegraphics[width=0.45\textwidth]{"
                + outplots[i + 2]
                + r"}"
                + "\n"
            )
            f.write(
                r"\includegraphics[width=0.45\textwidth]{"
                + outplots[i + 3]
                + r"}"
                + "\n"
            )
            f.write(r"\end{figure}" + "\n")
            f.write(r"\clearpage" + "\n")
        f.write(r"\end{document}")
    os.chdir(img.outdir)
    os.system("pdflatex -interaction=batchmode {0}".format(fname))
    os.chdir("..")
    img.logger.info("Done.")


def lineplot(img):
    """
    Generate PDF of channel cube diagnostic plots

    Inputs:
        img :: Imaging object
            The Imaging object

    Returns: Nothing
    """
    # Get center pixels
    center_x = int(img.cp["imsize"][0] / 2)
    center_y = int(img.cp["imsize"][1] / 2)

    # Loop over spws
    goodplots = []
    for spw in img.line_spws:
        # check that this spectral window exists
        fname = "{0}.spw{1}.{2}.channel.clean.image.fits".format(
            img.field, spw, img.stokes
        )
        if img.uvtaper:
            fname = "{0}.spw{1}.{2}.channel.uvtaper.clean.image.fits".format(
                img.field, spw, img.stokes
            )
        fname = os.path.join(img.outdir, fname)
        if not os.path.exists(fname):
            continue

        # Loop over all plot filenames
        if img.uvtaper:
            fitsfiles = [
                "{0}.spw{1}.{2}.channel.uvtaper.dirty.image.fits".format(
                    img.field, spw, img.stokes
                ),
                "{0}.spw{1}.{2}.channel.uvtaper.clean.image.fits".format(
                    img.field, spw, img.stokes
                ),
                "{0}.spw{1}.{2}.channel.uvtaper.clean.residual.fits".format(
                    img.field, spw, img.stokes
                ),
                "{0}.spw{1}.{2}.channel.uvtaper.clean.pbcor.image.fits".format(
                    img.field, spw, img.stokes
                ),
            ]
            maskfile = "{0}.spw{1}.{2}.channel.uvtaper.clean.mask.fits".format(
                img.field, spw, img.stokes
            )
            titles = [
                "{0} - {1} - Taper/Dirty".format(img.field, spw),
                "{0} - {1} - Taper/Clean".format(img.field, spw),
                "{0} - {1} - Taper/Residual".format(img.field, spw),
                "{0} - {1} - Taper/PBCorr".format(img.field, spw),
            ]
            labels = [
                "Flux Density (Jy/beam)",
                "Flux Density (Jy/beam)",
                "Flux Density (Jy/beam)",
                "Flux Density (Jy/beam)",
            ]
            vlims = [
                (None, None),
                (None, None),
                (None, None),
                (None, None),
            ]
        else:
            fitsfiles = [
                "{0}.spw{1}.{2}.channel.dirty.image.fits".format(
                    img.field, spw, img.stokes
                ),
                "{0}.spw{1}.{2}.channel.clean.image.fits".format(
                    img.field, spw, img.stokes
                ),
                "{0}.spw{1}.{2}.channel.clean.residual.fits".format(
                    img.field, spw, img.stokes
                ),
                "{0}.spw{1}.{2}.channel.clean.pbcor.image.fits".format(
                    img.field, spw, img.stokes
                ),
            ]
            maskfile = "{0}.spw{1}.{2}.channel.clean.mask.fits".format(
                img.field, spw, img.stokes
            )
            titles = [
                "{0} - {1} - Dirty".format(img.field, spw),
                "{0} - {1} - Clean".format(img.field, spw),
                "{0} - {1} - Residual".format(img.field, spw),
                "{0} - {1} - PBCorr".format(img.field, spw),
            ]
            labels = [
                "Flux Density (Jy/beam)",
                "Flux Density (Jy/beam)",
                "Flux Density (Jy/beam)",
                "Flux Density (Jy/beam)",
            ]
            vlims = [
                (None, None),
                (None, None),
                (None, None),
                (None, None),
            ]

        # Generate figures for each Stokes parameter
        for ind, stokes in enumerate(img.stokes):
            # Loop over figures
            for fitsfile, title, label, vlim in zip(fitsfiles, titles, labels, vlims):

                # Open fits file, generate WCS
                hdulist = fits.open(os.path.join(img.outdir, fitsfile))
                hdu = hdulist[0]
                wcs = WCS(hdu.header)

                # Generate figure
                plt.ioff()
                fig = plt.figure()
                ax = plt.subplot(projection=wcs.sub(["celestial"]))
                ax.set_title("{0} - {1}".format(title, stokes))
                center_chan = hdu.data.shape[1] / 2
                cax = ax.imshow(
                    hdu.data[ind, center_chan],
                    origin="lower",
                    interpolation="none",
                    cmap="binary",
                    vmin=vlim[0],
                    vmax=vlim[1],
                )
                # ax.grid(True,color='black',ls='solid')
                if img.cp["frame"] == "J2000":
                    ax.coords[0].set_major_formatter("hh:mm:ss")
                    ax.set_xlabel("RA (J2000)")
                    ax.set_ylabel("Declination (J2000)")
                elif img.cp["frame"] == "GALACTIC":
                    ax.set_xlabel("Galactic Longitude")
                    ax.set_ylabel("Galactic Latitude")

                # Plot clean mask as contour
                maskhdu = fits.open(os.path.join(img.outdir, maskfile))[0]
                _ = ax.contour(
                    maskhdu.data[ind, center_chan],
                    levels=[0.5],
                    origin="lower",
                    colors="r",
                    linewidths=2,
                )

                # Plot beam, if it is defined
                if "BMAJ" in hdu.header.keys():
                    cell = float(img.cp["cell"].replace("arcsec", ""))
                    beam_maj = hdu.header["BMAJ"] * 3600.0 / cell
                    beam_min = hdu.header["BMIN"] * 3600.0 / cell
                    beam_pa = hdu.header["BPA"]
                    ellipse = Ellipse(
                        (
                            center_x - int(3.0 * center_x / 4),
                            center_y - int(3.0 * center_y / 4),
                        ),
                        beam_min,
                        beam_maj,
                        angle=beam_pa,
                        fill=True,
                        zorder=10,
                        hatch="///",
                        edgecolor="black",
                        facecolor="white",
                    )
                    ax.add_patch(ellipse)
                elif len(hdulist) > 1:
                    hdu = hdulist[1]
                    cell = float(img.cp["cell"].replace("arcsec", ""))
                    beam_maj = hdu.data["BMAJ"][center_chan] / cell
                    beam_min = hdu.data["BMIN"][center_chan] / cell
                    beam_pa = hdu.data["BPA"][center_chan]
                    ellipse = Ellipse(
                        (
                            center_x - int(3.0 * center_x / 4),
                            center_y - int(3.0 * center_y / 4),
                        ),
                        beam_min,
                        beam_maj,
                        angle=beam_pa,
                        fill=True,
                        zorder=10,
                        hatch="///",
                        edgecolor="black",
                        facecolor="white",
                    )
                    ax.add_patch(ellipse)

                # Plot colorbar
                cbar = fig.colorbar(cax, fraction=0.046, pad=0.04)
                cbar.set_label(label)

                # Re-scale to fit, then save
                fname = fitsfile.replace(".fits", ".{0}.pdf".format(stokes))
                fig.savefig(os.path.join(img.outdir, fname), bbox_inches="tight")
                plt.close(fig)
                plt.ion()
                goodplots.append(fname)

            # Generate spectrum
            if img.uvtaper:
                fitsfile = "{0}.spw{1}.{2}.channel.uvtaper.clean.image.fits".format(
                    img.field, spw, img.stokes
                )
            else:
                fitsfile = "{0}.spw{1}.{2}.channel.clean.image.fits".format(
                    img.field, spw, img.stokes
                )
            hdu = fits.open(os.path.join(img.outdir, fitsfile))[0]
            spec = hdu.data[ind, :, center_x, center_y]
            isnan = spec == 0.0
            spec[isnan] = np.nan
            velo = (
                np.arange(len(spec)) * hdu.header["CDELT3"] + hdu.header["CRVAL3"]
            ) / 1000.0

            # Generate figure
            plt.ioff()
            fig = plt.figure()
            ax = plt.subplot()
            ax.plot(velo, spec, "k-")
            ax.set_xlabel("Velocity (km/s)")
            ax.set_ylabel("Flux Density (Jy/beam)")
            ax.set_xlim(np.nanmin(velo), np.nanmax(velo))
            ybuff = 0.1 * (np.nanmax(spec) - np.nanmin(spec))
            ax.set_ylim(np.nanmin(spec) - ybuff, np.nanmax(spec) + ybuff)
            if img.uvtaper:
                ax.set_title(
                    "{0} - {1} - Taper/Center - {2}".format(img.field, spw, stokes)
                )
            else:
                ax.set_title("{0} - {1} - Center - {2}".format(img.field, spw, stokes))
            ax.grid(False)
            fig.tight_layout()
            fname = fitsfile.replace(".fits", ".{0}.spec.pdf".format(stokes))
            fig.savefig(os.path.join(img.outdir, fname), bbox_inches="tight")
            plt.close(fig)
            plt.ion()
            goodplots.append(fname)

    # Generate PDF of plots
    # need to fix filenames so LaTeX doesn't complain
    goodplots = ["{" + fn.replace(".pdf", "") + "}.pdf" for fn in goodplots]
    img.logger.info("Generating PDF...")
    fname = "{0}.lineplots.tex".format(img.field)
    if img.uvtaper:
        fname = "{0}.uvtaper.lineplots.tex".format(img.field)
    with open(os.path.join(img.outdir, fname), "w") as f:
        f.write(r"\documentclass{article}" + "\n")
        f.write(r"\usepackage{graphicx}" + "\n")
        f.write(r"\usepackage[margin=0.1cm]{geometry}" + "\n")
        f.write(r"\begin{document}" + "\n")
        for i in range(0, len(goodplots), 5):
            f.write(r"\begin{figure}" + "\n")
            f.write(r"\centering" + "\n")
            f.write(r"\includegraphics[width=0.45\textwidth]{" + goodplots[i] + "}\n")
            f.write(
                r"\includegraphics[width=0.45\textwidth]{"
                + goodplots[i + 1]
                + r"} \\"
                + "\n"
            )
            f.write(
                r"\includegraphics[width=0.45\textwidth]{" + goodplots[i + 2] + "}\n"
            )
            f.write(
                r"\includegraphics[width=0.45\textwidth]{"
                + goodplots[i + 3]
                + r"} \\"
                + "\n"
            )
            f.write(
                r"\includegraphics[width=0.45\textwidth]{" + goodplots[i + 4] + "}\n"
            )
            f.write(r"\end{figure}" + "\n")
            f.write(r"\clearpage" + "\n")
        f.write(r"\end{document}")
    os.chdir(img.outdir)
    os.system("pdflatex -interaction=batchmode {0}".format(fname))
    os.chdir("..")
    img.logger.info("Done.")
