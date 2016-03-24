#! /bin/tcsh

set currdir = `pwd`

    echo "running TSF SmallSet"
    ${currdir}/TryTwoStateFit.py All SmallSet 
    echo "running OSF SmallSet"
    ${currdir}/TryOneStateFit.py All SmallSet 
    echo "running Summation"
    ${currdir}/TrySummation.py SmallSet
    
