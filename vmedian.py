"""
Vector Median
"""
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2023 Quiet Skies LLC
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
# HISTORY
# 19JUN20 GIL try to optimize 
#

import numpy as np

def vmedian( ins, n2):
    """
    Vector Median takes an array of float values and returns the median of 
    of that vector, with half width n2
    """

    n = len( ins)
    outs = ins             # initialize the output values
    
    for i in range(n):
        i0 = max( 0, i - n2)
        i1 = min( i + n2 + 1, n)
        outs[i] = np.median( ins[i0:i1])
    return outs
    # end vmedian()

def vsmooth( ins, n2):
    """
    Vector Median takes an array of float values and returns the median of 
    of that vector, with half width n2
    """

    n = len( ins)
    outs = ins             # initialize the output values
    
    for i in range(n):
        i0 = max( 0, i - n2)
        i1 = min( i + n2, n)
        outs[i] = np.average( ins[i0:i1])
    return outs
    # end vmedian()
    
