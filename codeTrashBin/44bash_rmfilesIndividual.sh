#!/bin/bash


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Note: the < #   PAIR_<ndx> > is used to identify training blocks
#   version 0.10  :  cleaning done folder
# 
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


# directory to be processed, if it is void, do all & final cleaning
FOLDER_INDEX=1

# Files to be used as reference
MAEFILE=MAE_PAIR_total.txt
MAEFILTER=0.12

# search keyword in MAEFILTER, finally in a combo:
# inside MAEFILTER:   ${MAE_KEYWORD}${FOLDER_INDEX},      e.g.: ChargePair_15-6
# for folder      :   ${FOLDER_KEYWORD}${FOLDER_INDEX},   e.g.: Training_Topfile_15-6

MAE_KEYWORD='ChargePair_'
FOLDER_KEYWORD='Training_Topfile_'


#
# END of user inputs
#

FOLDER=${FOLDER_KEYWORD}${FOLDER_INDEX}
SEARCH_KEYWORD=${MAE_KEYWORD}${FOLDER_INDEX}

if [[ ! -f $MAEFILE ]]; then { echo "Error: MAEFILE < $MAEFILE > is not found"; exit 1; } fi
if [[ ! -d $FOLDER ]]; then { echo "Error: FOLDER < $FOLDER > is not found"; exit 1; } fi
if (( $(ls $FOLDER | wc -l) == 0 )); then { echo "Error: FOLDER < $FOLDER > is empty"; exit 1; } fi

# Now get line-nm
ndxbeg=$(sed -n "/\#\ *${SEARCH_KEYWORD}/=" $MAEFILE)
if [[ -z "$ndxbeg" ]]; then { echo "Error: cannot get line-nm in MAEFILE"; exit 1; } fi


cwd=$(pwd)
mol=${cwd//*\/}
echo "Note: Processing molecule: < $mol >"
echo "Note: Do you want to continue? y/yes, else not"

read tmp
if [[ -z "$tmp" || ! ( "$tmp" == 'yes' || "$tmp" == 'y' ) ]]
then
    echo "You decided to quit..."
    exit 0
fi

if [[ -z $(ls | grep '.err') ]]
then
    echo "Warning: it seems the simulation is NOT done yet."
    echo "Do you REALLY want to continue? y/yes, else not"
    read tmp
    if [[ -z "$tmp" || ! ( "$tmp" == 'yes' || "$tmp" == 'y' ) ]]
    then
        echo "You decided to quit..."
        exit 0
    fi
fi


# Now get each Training PAIR numbers
ndxend=$(sed -n "$ndxbeg,\$p" $MAEFILE | sed -n "/\#\ +${MAE_KEYWORD}/=" | head -n 1)
if [[ -z "$ndxend" ]]; then ndxend=$(cat $MAEFILE | wc -l); fi


# lenndx: the length of each PAIR
content=($(sed -n "${ndxbeg},${ndxend}p" $MAEFILE | grep '^PAIR' | head -n 1))
lenndx=${#content[*]}

if (( $lenndx == 0 )); then { echo "Fatal Error: invalid input file < $MAEFILE >"; exit 1; } fi


KEYWORD=MAE
echo "Note: Setting processing file to < $MAEFILE >..."
echo "Note: Setting Keyword to < $KEYWORD >..."
echo ''

# ndx: the index number of MAE value in echo PAIR line
for ((ndx=0; $ndx < $lenndx; ndx++))
do
    if [[ ${content[$ndx]} == $KEYWORD ]]; then break; fi
done
((ndx++))

if (($ndx >= $lenndx)); then { echo "Fatal Error: No KEYWORD=$KEYWORD is found"; exit 1; } fi


data=''
content=($(sed -n "${ndxbeg},${ndxend}p" $MAEFILE | grep '^PAIR'))
lentot=${#content[*]}
for ((i=$ndx; $i<$lentot; i=$i+$lenndx))
do
    tmp=$(echo ${content[$i]} | tr [A-Z] [a-z])
    if [[ $tmp != 'nan' && ! $tmp =~ ^[+-]?([0-9]*[.])?([0-9]+)?$ ]]
    then
        echo "Fatal Error: Invalid $KEYWORD value in line"
        for ((j=$i-$ndx; $j < $i - $ndx + $lenndx; j++)); do echo -n "${content[$j]}  "; done
        echo ''
        exit 1
    fi
    data="$data $tmp"
done


# double check folder, warning on the multi-processing
tmp=($data)
cnt=0
for i in $(ls $FOLDER); do { if [[ -d $i ]]; then ((cnt++)); fi } done
if (( ${#tmp[*]} > $cnt ))
then
    echo "Warning: FOLDER < $FOLDER > seems to have been processed, do you want to continue?"
    read ts
    if [[ -z "$ts" || ! ( "$ts" == 'yes' || "$ts" == 'y' ) ]]
    then
        echo "You decided to quit..."
        exit 0
    fi
fi


# Process data: only MAE less than MAEFILTER will be left
# refstr  :  index number to be removed
# Tip: let < cnt > starting at zero
cnt=0
refstr=''
for i in $data
do
    ((cnt++))
    if [[ $i == 'nan' ]]; then { refstr="$refstr $cnt"; continue; } fi

    tmp=$(echo "$MAEFILTER - $i" | bc -l | grep '-')
    if [[ -n "$tmp" ]]; then refstr="$refstr $cnt"; fi
done

cd $FOLDER
for i in $refstr
do
    if [[ -d dir_$i ]]
    then
        echo "Note: Deep Processing dir_$i"
        rm -rf dir_$i
    fi
done

for i in dir_*
do
    if [[ -d $i ]]
    then
        cd $i
        echo "Note: Processing < $i >"
        rm -f energy* \#energy* mdout* \#mdout* step* \#step*
        rm -f min* npt*
        rm -f prod.cpt prod.log prod_surftmp.gro
        rm -f surf.cpt surf.log tmp*

        for j in vis_*
        do
            if [[ -d $j ]]
            then
                cd $j
                rm -f energy* \#energy* mdout* \#mdout* step* \#step*
                rm -f vis.cpt vis.log vis.trr
                cd ../
            fi
        done

        cd ../
    fi
done



echo ''
echo "Everything is DONE"
