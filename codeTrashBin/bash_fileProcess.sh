#!/bin/bash

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# version 0.20  :  process files, mkdir; and rename them to MOL.*
# version 0.30  :  backup original file and prepare for simulation
#

# contents in FILE to be processed
FILE=files


# Folder storing original files
FLPG=LigParGen

# Folder storing processed files
FMOL=FMOL

# run files in copy folders
FCOPY=SelectedFiles


#
# End of user inputs
#

nmargv=$#

FILE=${FILE:=files}
if (( $nmargv == 0 ))
then
    :
elif (( $nmargv == 1 ))
then
    if [[ -f $1 ]]
    then
        FILE=$1
    else
        echo "Error: input is not a file"
        exit 1
    fi
else
    echo "Error: wrong number of parameters"
    exit 1
fi

if [[ ! -f $FILE ]]; then { echo "Note: FILE < $FILE > does not exist"; exit 1; } fi
if [[ ! -d $FMOL ]]; then mkdir $FMOL; fi
if [[ ! -d $FLPG ]]; then mkdir $FLPG; fi
if (( $(ls $FMOL | wc -l) > 0 )); then { echo "Folder < $FMOL > has to be empty"; exit 1; } fi


content=''
for f in $(cat $FILE)
do
    if [[ -f $f ]]
    then
        content="$content $f"
    else
        echo "Warning: file < $f > does not exist"
    fi
done
content=($(echo "$content" | sort))
unicontent=""
for ((i=0; $i<${#content[*]}; i=$i+2))
do
    bo=true
    fa=${content[$i]}
    fa_name=${fa%.*}
    fa_ext=${fa/*.}

    ((j=$i+1))
    fb=${content[$j]}
    fb_name=${fb%.*}
    fb_ext=${fb/*.}

    if [[ "$fa_name" != "$fb_name" ]]; then { echo "Warning: file errors for < $fa > and < $fb >"; continue; } fi
    if [[ ( "$fa_ext" == 'itp' || "$fb_ext" == 'pdb' ) || ( "$fa_ext" == "pdb" || "$fb_ext" == 'itp' ) ]]
    then
        unicontent="$unicontent $fa_name"
    else
        echo "Warning: file errors for < $fa > and < $fb >"
    fi
done


# backup and processing
for f in $unicontent
do
    itp=${f}.itp
    pdb=${f}.pdb

    echo "Note: Processing < ${f} itp/pdb >"

    if [[ -d $FLPG/$f ]]
    then
        echo "Warning: Folder < $FLPG/$f > already exists"
        cnt=1
        while true
        do
            if [[ ! -d $FLPG/${f}_$cnt ]]; then break; fi
            ((cnt++))
        done
        foriginal=$FLPG/${f}_$cnt
    else
        foriginal=$FLPG/${f}
    fi
    echo "Note: Moving files < $f > to Folder < $foriginal >"
    mkdir $foriginal
    cp $itp $pdb $foriginal

    sed -i 's/UNK/MOL/' $itp
    sed -i 's/UNK/MOL/' $pdb

    sed -i '/^CONECT/d; /^TER/d' $pdb

    nm=($(sed -n '/\[ pairs \]/=' $itp))
    if [[ -n "${nm[0]}" ]]
    then
        sed -i "${nm[0]},\$d" $itp
    fi


    mkdir $FMOL/$f
    mv $itp $FMOL/$f/MOL.itp
    mv $pdb $FMOL/$f/MOL.pdb

    cp $FCOPY/* $FMOL/$f

    echo "$f"  >> NEWfiles
done

echo '' >> NEWfiles
echo '' >> NEWfiles
echo '' >> NEWfiles

echo ''
echo "Done Everything"
