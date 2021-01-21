"""
utils.py - WISP Utilities

General utilities.

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
Trey V. Wenger August 2019 - V2.1
    Improve code readability.
"""

import os
import re
import glob
import shutil
import numpy as np

import casac
from casac import *
from tasks import imregrid
from numpy import fabs


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
    def convert(text):
        if text.isdigit():
            return int(text)
        return text.lower()

    # define the sorting algorithm
    def alphanum_key(key):
        return [convert(c) for c in re.split("([0-9]+)", key)]

    # return the sorted list
    return sorted(mylist, key=alphanum_key)


def generate_pdf(location):
    """
    Compile images into a single PDF.

    Inputs: location
      location :: string
        The directory containing the images. Also defines the
        PDF filename (`location`.pdf)

    Returns: Nothing
    """
    iplot = 0

    # Get images
    fnames = glob.glob("{0}/*.png".format(location))
    fnames = natural_sort(fnames)

    # Generate TeX document
    pdfname = "{0}.tex".format(location)
    with open(pdfname, "w") as fout:
        fout.write(r"\documentclass{article}" + "\n")
        fout.write(r"\usepackage{graphicx,subfig}" + "\n")
        fout.write(r"\usepackage[margin=0.1cm]{geometry}" + "\n")
        fout.write(r"\begin{document}" + "\n")
        fout.write(r"\begin{figure}" + "\n")
        fout.write(r"\centering" + "\n")
        for fname in fnames:
            if iplot > 0 and iplot % 6 == 0:
                fout.write(r"\end{figure}" + "\n")
                fout.write(r"\clearpage" + "\n")
                fout.write(r"\begin{figure}" + "\n")
                fout.write(r"\centering" + "\n")
            elif iplot > 0 and iplot % 2 == 0:
                fout.write(r"\end{figure}" + "\n")
                fout.write(r"\begin{figure}" + "\n")
                fout.write(r"\centering" + "\n")
            fout.write(
                r"\includegraphics[width=0.45\textwidth]{" + fname + "}" + "\n"
            )
            iplot += 1
        fout.write(r"\end{figure}" + "\n")
        fout.write(r"\end{document}" + "\n")

    # Compile PDF
    os.system("pdflatex -interaction=batchmode {0}.tex".format(location))


def get_smodels(cal):
    """
    Get secondary calibrator polarization estimates from calibrated data.

    Inputs:
        cal :: Calibration object
            The calibration object

    Returns: smodels
        smodels :: Dictionary
            Keys are secondary calibrator field names, values are
            [1.0, Q/I, U/I, 0.0] estimated from corrected_data
    """
    smodels = {}

    # loop over all calibrators, but only save secondary calibrators
    for field in cal.calibrators:
        # get fieldid(s)
        cal.casa.msmd.open(cal.vis)
        fieldids = list(cal.casa.msmd.fieldsforname(field))
        cal.casa.msmd.close()

        # estimate Q and U for each continuum spw
        # average to get smodel
        spw_q = []
        spw_u = []
        for spw in cal.cross_spws:
            # get data
            cal.casa.ms.open(cal.vis)
            cal.casa.ms.selectinit(datadescid=int(spw))
            cal.casa.ms.select({"field_id": fieldids})
            cal.casa.ms.selectpolarization(["I", "Q", "U", "V"])
            data = cal.casa.ms.getdata(["corrected_data", "flag"])
            cal.casa.ms.close()

            # get un-flagged data
            i = data["corrected_data"][0, :, :][~data["flag"][0, :, :]]
            q = data["corrected_data"][1, :, :][~data["flag"][1, :, :]]
            u = data["corrected_data"][2, :, :][~data["flag"][2, :, :]]
            v = data["corrected_data"][3, :, :][~data["flag"][3, :, :]]
            ave_i = np.complex(np.mean(np.real(i)), np.mean(np.imag(i)))
            ave_q = np.complex(np.mean(np.real(q)), np.mean(np.imag(q)))
            ave_u = np.complex(np.mean(np.real(u)), np.mean(np.imag(u)))
            ave_v = np.complex(np.mean(np.real(v)), np.mean(np.imag(v)))
            abs_i = np.sign(np.real(ave_i)) * np.abs(ave_i)
            abs_q = np.sign(np.real(ave_q)) * np.abs(ave_q)
            abs_u = np.sign(np.real(ave_u)) * np.abs(ave_u)
            abs_v = np.sign(np.real(ave_v)) * np.abs(ave_v)
            cal.logger.info(
                "Field %s spw %s has (I,Q,U,V) = (%.4f, %.4f, %.4f, %.4f)",
                field,
                spw,
                abs_i,
                abs_q,
                abs_u,
                abs_v,
            )
            spw_q.append(abs_q / abs_i)
            spw_u.append(abs_u / abs_i)

        # average spws
        estimate_q = np.nanmean(spw_q)
        std_q = np.nanstd(spw_q)
        estimate_u = np.nanmean(spw_u)
        std_u = np.nanstd(spw_u)
        cal.logger.info(
            "Field %s spw average: Q/I=%.4f U/I=%.4f",
            field,
            estimate_q,
            estimate_u,
        )
        cal.logger.info(
            "Field %s spw std. dev.: Q/I=%.4f U/I=%.4f",
            field,
            std_q,
            std_u,
        )

        # skip if this is a polarization and/or flux calibrator
        # and therefore already has a model from setjy
        if field in cal.flux_cals + cal.pol_leak_cals + cal.pol_angle_cals:
            cal.logger.info(
                "Not saving smodel because this is a "
                "flux or polarization calibrator."
            )
            continue
        smodels[field] = [1, estimate_q, estimate_u, 0]
    return smodels


def makePB(
    vis="",
    field="",
    spw="",
    timerange="",
    uvrange="",
    antenna="",
    observation="",
    intent="",
    scan="",
    imtemplate="",
    outimage="",
    pblimit=0.2,
    stokes="I",
):
    """
    (Modified from recipes.makepb, added multiple stokes capability)

    Make a PB image using the imager tool, onto a specified image coordinate
    system

    This function can be used along with tclean to make .pb images for
    gridders that do not already do it (i.e. other than mosaic, awproject)

    This script takes an image to use as a template coordinate system,
    attempts to set up an identical coordinate system with the old imager tool,
    makes a PB for the telescope listed in the MS observation subtable, and
    regrids it (just in case) to the target coordinate system). This can be
    used for single fields and mosaics.
    """
    tb = casac.table()
    im = casac.imager()
    ia = casac.image()
    me = casac.measures()
    qa = casac.quanta()
    tb.open(vis + "/OBSERVATION")
    tel = tb.getcol("TELESCOPE_NAME")[0]
    tb.close()
    tb.open(vis + "/SPECTRAL_WINDOW")
    mfreqref = tb.getcol("MEAS_FREQ_REF")[0]
    tb.close()
    if mfreqref == 64:
        print(
            "MAKEPB : This function is using old imager tool, "
            "Undefined frame may not be handled properly."
        )
    ia.open(imtemplate)
    csysa = ia.coordsys()
    csys = csysa.torecord()
    shp = ia.shape()
    ia.close()
    dirs = csys["direction0"]
    phasecenter = me.direction(
        dirs["system"],
        qa.quantity(dirs["crval"][0], dirs["units"][0]),
        qa.quantity(dirs["crval"][1], dirs["units"][1]),
    )
    cellx = qa.quantity(fabs(dirs["cdelt"][0]), dirs["units"][0])
    celly = qa.quantity(fabs(dirs["cdelt"][1]), dirs["units"][1])
    nchan = shp[3]
    # assumes refpix is zero
    start = qa.quantity(csysa.referencevalue()["numeric"][3], csysa.units()[3])
    mestart = me.frequency("LSRK", start)
    step = qa.quantity(csysa.increment()["numeric"][3], csysa.units()[3])
    smode = "mfs"
    if nchan > 1:
        smode = "frequency"
    im.open(vis)
    im.selectvis(
        field=field,
        spw=spw,
        time=timerange,
        intent=intent,
        scan=scan,
        uvrange=uvrange,
        baseline=antenna,
        observation=observation,
    )
    im.defineimage(
        nx=shp[0],
        ny=shp[0],
        phasecenter=phasecenter,
        cellx=qa.tos(cellx),
        celly=qa.tos(celly),
        stokes=stokes,
        nchan=nchan,
        start=mestart,
        step=step,
        mode=smode,
    )
    im.setvp(dovp=True, telescope=tel)
    im.makeimage(type="pb", image=outimage + ".tmp")
    im.close()
    if mfreqref == 64:  # skip this step if the frame is 'Undefined'
        shutil.copytree(outimage + ".tmp", outimage)
    else:
        print("MAKEPB : Regrid to desired coordinate system")
        imregrid(
            imagename=outimage + ".tmp",
            template=imtemplate,
            output=outimage,
            overwrite=True,
            asvelocity=False,
        )
    shutil.rmtree(outimage + ".tmp")
    ia.open(outimage)
    ia.calcmask("'" + outimage + "'>" + str(pblimit))
    ia.close()
