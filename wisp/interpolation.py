"""
interpolation.py - Interpolate through some channels in a measurement set

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

import numpy as np
from scipy.interpolate import interp1d

from .utils import natural_sort


def interpolate_channels(cal, spws, chans):
    """
    Edit the measurement set to replace bad channels with
    interpolated values. Simple linear interpolation between
    neighbors, in both phase and amplitude.

    Inputs:
        cal :: Calibration object
            The calibration object
        spws :: list of integers
            The spectral windows to use
        chans :: list of integers
            The channels through which to interpolate

    Returns: Nothing
    """
    bad_chans = np.array(chans)

    # Open MS for modifications, and get list of data_desc_ids
    cal.casa.ms.open(cal.vis, nomodify=False)
    spwinfo = cal.casa.ms.getspectralwindowinfo()
    datadescids = natural_sort(spwinfo.keys())
    for datadescid in datadescids:
        # Check that the spw associated with this data_desc_id
        # is one that needs interpolated
        if spwinfo[datadescid]["SpectralWindowId"] not in spws:
            continue
        cal.logger.info("Working on spw %d", spwinfo[datadescid]["SpectralWindowId"])
        nchans = spwinfo[datadescid]["NumChan"]
        chans = np.arange(nchans)
        mask = np.zeros(nchans, dtype=bool)
        mask[bad_chans] = True

        # Select data_desc_id, and initialize iterator
        cal.casa.ms.selectinit(datadescid=int(datadescid))
        cal.casa.ms.iterinit()
        cal.casa.ms.iterorigin()

        # Iterate over chunks
        while True:

            # get data
            rec = cal.casa.ms.getdata(["data"])

            # interpolate real part
            real = np.real(rec["data"])
            real_interp = interp1d(chans[~mask], real[:, ~mask, :], axis=1)
            real[:, mask, :] = real_interp(chans[mask])

            # interpolate phase
            imag = np.imag(rec["data"])
            imag_interp = interp1d(chans[~mask], imag[:, ~mask, :], axis=1)
            imag[:, mask, :] = imag_interp(chans[mask])

            # save and store
            rec["data"] = real + 1.0j * imag
            cal.casa.ms.putdata(rec)

            # Get next chunk
            if not cal.casa.ms.iternext():
                break

        # Terminate iterator and reset selection
        cal.casa.ms.iterend()
        cal.casa.ms.selectinit(reset=True)
        cal.logger.info("Done.")

    # Close measurement set
    cal.casa.ms.close()
