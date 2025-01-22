#!/bin/bash
# 
# Simple script to read .dot file, find out whether it encodes a directed or
# undirected graph, and launch the corresponding GraphViz utility to generate
# images in png and svg format.
# 
# 
# Copyright (c) 2022-2024 Vasilios Raptis <v.raptis@external.euc.ac.cy>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# 
#-------------------------------------------------------------------------------
# 
dotfil=$1
svgfil=`basename $dotfil dot`"svg"
pngfil=`basename $dotfil dot`"png"
echo $svgfil $pngfil
#dotcmd="neato -Goverlap=scale -Gsplines=true -Gstart=$RANDOM"
#dotcmd="neato -Goverlap=scale -Gsplines=true -Gstart=3"   # this seed produced nicer schemes
f=`grep neato $1`
echo $f
if [[ $f = '' ]]
then
    dotcmd="neato -Goverlap=scale -Gsplines=true"
else
    dotcmd="dot   -Goverlap=scale -Gsplines=true"
fi
eval $dotcmd -Tpng $dotfil > $pngfil
eval $dotcmd -Tsvg $dotfil > $svgfil
eog $pngfil
