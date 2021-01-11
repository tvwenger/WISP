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

__version__ = '2.1'

import os
import re
import glob


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
        return [convert(c) for c in re.split('([0-9]+)', key)]

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
    #
    # Get images
    #
    fnames = glob.glob('{0}/*.png'.format(location))
    fnames = natural_sort(fnames)
    #
    # Generate TeX document
    #
    pdfname = '{0}.tex'.format(location)
    with open(pdfname, 'w') as fout:
        fout.write(r'\documentclass{article}'+'\n')
        fout.write(r'\usepackage{graphicx,subfig}'+'\n')
        fout.write(r'\usepackage[margin=0.1cm]{geometry}'+'\n')
        fout.write(r'\begin{document}'+'\n')
        fout.write(r'\begin{figure}'+'\n')
        fout.write(r'\centering'+'\n')
        for fname in fnames:
            if iplot > 0 and iplot % 6 == 0:
                fout.write(r'\end{figure}'+'\n')
                fout.write(r'\clearpage'+'\n')
                fout.write(r'\begin{figure}'+'\n')
                fout.write(r'\centering'+'\n')
            elif iplot > 0 and iplot % 2 == 0:
                fout.write(r'\end{figure}'+'\n')
                fout.write(r'\begin{figure}'+'\n')
                fout.write(r'\centering'+'\n')
            fout.write(
                r'\includegraphics[width=0.45\textwidth]{'+fname+'}'+'\n')
            iplot += 1
        fout.write(r'\end{figure}'+'\n')
        fout.write(r'\end{document}'+'\n')
    #
    # Compile PDF
    #
    os.system('pdflatex -interaction=batchmode {0}.tex'.format(location))
