#!/bin/bash

# copy bbl file to the previous version
# cp latex/CIJ_JGR.bbl ../RadiationPatternAtlas3/latex/

latexdiff --flatten \
--append-context2cmd="twplot,fplot,plot,tplot,ddplot" \
--append-mboxsafecmd="cite,citep" \
--graphics-markup="1" \
--allow-spaces \
 --disable-auto-mbox \
 --add-to-config "MATHENV=beq;eeq" \
../RadiationPatternAtlas3/latex/CIJ_JGR.tex latex/CIJ_JGR.tex > latex/diffW3.tex

# --exclude-safecmd="beq,eeq" \
# --show-config \
#--show-config \
#--exclude-safecmd="twplot,fplot,plot,tplot,ddplot,dplot" \
#--exclude-mboxsafecmd="twplot,fplot,plot,tplot,ddplot,dplot" \
#--append-textcmd="twplot,fplot,plot,tplot,ddplot,dplot" \
#--append-mboxsafecmd="cite,citep" \
#--exclude-safecmd="citeauthoryear,bibitem" \