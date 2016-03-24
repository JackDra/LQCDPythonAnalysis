#! /bin/tcsh

set currdir = `pwd`

set currentlist = 'Scalar PsVector giDi Vector'
# set currentlist = 'Vector PsVector giDi'

echo "running TwoPt"
${currdir}/RunCorrAnaCurr.py Scalar TwoPt

foreach icurr ( ${currentlist} )
    
    echo "running ${icurr} CM 29"
    ${currdir}/RunCorrAnaCurr.py ${icurr} CM 29 
    # echo "running ${icurr} REvec 26"
    # ${currdir}/RunCorrAnaCurr.py ${icurr} REvec 26
    # echo "running ${icurr} TSink '26 32 35 38'"
    # ${currdir}/RunCorrAnaCurr.py ${icurr} TSink '26 32 35 38'

end
