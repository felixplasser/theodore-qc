#!/bin/bash

PDIR=`pwd`

echo "theo_test.bash [<module>]"
echo "  available modules: standard, all, cclib"

echo "Starting theo_test.bash"
echo "THEODIR=$THEODIR"

stddirs="pyrrole.qcadc hexatriene.colmrci fa2.ricc2 pv2p.escf pv2p.qctddft ir_c3n3.qctddft pyridine.ricc2"
cclibdirs="fa2.cclib"

if [ $# == 0 ]
then
    rundirs=$stddirs
else
    case "$1" in
    "all") rundirs="$stddirs $cclibdirs";;
    "standard") rundirs="$stddirs";;
    "cclib") rundirs="$cclibdirs";;    
    esac
fi

rm -r RUN_THEO_TEST
mkdir RUN_THEO_TEST || exit 1

tchk=0
for dir in $rundirs
do
    echo
    echo "================================================"
    echo
    echo "Starting test $dir ..."
    sdir="$THEODIR/EXAMPLES/$dir"
    
    rdir="$PDIR/RUN_THEO_TEST/$dir"
    if [ -d $rdir ]
    then
        echo " ERROR:"
        echo "$rdir already exists! Please delete it or run in a different directory."
        exit 5
    fi
    
    cp -r "$sdir/QC_FILES" $rdir
    cd $rdir
    
    chk=0
    for ifile in `ls "$sdir/IN_FILES"`
    do
        cp $sdir/IN_FILES/$ifile .
        
        dtype=`echo $ifile | cut -d '.' -f 1`
        atype=`echo $ifile | cut -d '.' -f 3`
        comm="analyze_$dtype.py -ifile $ifile"
        echo $comm
        $comm > analyze_$atype.out
        lchk=$?
        if [ "$lchk" -ne 0 ]
        then
            echo "  ... failed!"
        fi
        
        chk=$((chk+lchk))
    done
    
    echo
    echo "Checking primary output files:"
    for rfile in `ls "$sdir/REF_FILES"`
    do
        echo "  -> $rfile"
        diff -I TheoDORE "$sdir/REF_FILES/$rfile" $rfile
        chk=$((chk+$?))
    done
    
    echo
    if [ -d $sdir/REF_FILES_SEC ]
    then
        echo "Checking secondary output files:"
        echo "  (These may show some numerical inaccuracies)"
        for rfile in `ls "$sdir/REF_FILES_SEC"`
        do
            echo "  -> $rfile"
            diff -q "$sdir/REF_FILES_SEC/$rfile" $rfile
            diff "$sdir/REF_FILES_SEC/$rfile" $rfile >> $PDIR/RUN_THEO_TEST/diff_sec.out
        done
    fi

    echo    
    echo " *** Test $dir finished (error code: $chk)."
    tchk=$((tchk+chk))
done
echo
echo "================================================"
echo
echo " *** All tests finished (number of errors: $tchk)"
echo
echo "================================================"

exit $tchk