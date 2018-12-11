#!/bin/bash

latexdiff --flatten --append-textcmd="twplot, fplot, tplot, ddplot, dplot" \
../RadiationPatternAtlas3/latex/CIJ_JGR.tex latex/CIJ_JGR.tex > latex/diffW3.tex

