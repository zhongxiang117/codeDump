#!/bin/bash


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Note: the < #   PAIR_ > is used to identify training blocks
#
#   Version 0.2   :   Deep cleaning based on MAE
#   Version 0.3   :   More robust dealing with resuming, check Help Usage
#   Versino 0.31  :   Suppress non-necessary printing in normal processing
#   Version 0.4   :   Remove empty folders, and reset MAE=0.2
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


# Files to be used as reference
FILE=MAE_PAIR_total.txt

#
# END of user inputs
#

FILE=${FILE:=MAE_PAIR_total.txt}

HELP_argument="Usage:    $0 [FILE]

>>> FILE : the input reference file. Default is < MAE_PAIR_total.txt >

Note:
    For input file, the entry < #  PAIR_Charge_nm.txt > is used to 
    identify < gennm >, and the < nm > will be extracted and compared
    with real Training Number. If anything goes wrong, Exit.
"

nmargv=$#
if (( $nmargv == 0 ))
then
    :
elif (( $nmargv == 1 ))
then
    if [[ "$1" == 'help' || "$1" == '-h' || "$1" == '--help' || "$1" == '-help' ]]
    then
        echo "$HELP_argument"
        exit 0
    fi
else
    echo 'Error: Number of input parameters were wrong, quitting...'
    echo ''
    exit 1
fi

bool_error=false
if (( $nmargv == 0 ))
then
    if ! [[ -f $FILE ]]; then bool_error=true; fi
else
    FILE=$1
    if ! [[ -f $1 ]]; then bool_error=true; fi
fi

if $bool_error
then
    echo "Error: < $FILE > file does not exist, quitting..."
    echo ''
    exit 1
fi


cwd=$(pwd)
mol=${cwd//*\/}
echo "Note: Processing molecule: $mol"
echo "Note: do you want to continue? y/yes, else not"

read tmp
if [[ -z "$tmp" || ! ( "$tmp" == 'yes' || "$tmp" == 'y' ) ]]
then
    echo "You decided to quit..."
    exit 0
fi

if [[ -z $(ls | grep '.err') ]]
then
    echo "Warning: it seems the simulation is not done"
    echo "Do you REALLY want to continue? y/yes, else not"
    read tmp
    if [[ -z "$tmp" || ! ( "$tmp" == 'yes' || "$tmp" == 'y' ) ]]
    then
        echo "You decided to quit..."
        exit 0
    fi
fi


KEYWORD=MAE
echo "Note: Setting processing file to < $FILE >..."
echo "Note: Setting Keyword to < $KEYWORD >..."



# Now get each Training PAIR numbers
ndxstr=$(sed -n '/\#\ *PAIR_/=' $FILE)
if [[ -z "$ndxstr" ]]; then { echo "Fatal Error: cannot get < gennm >"; exit 1; } fi

# Add the total number of lines
tot=$(cat $FILE | wc -l)
ndxlist=( $ndxstr $tot )


# lenndx: the length of each PAIR
content=($(sed -n "${ndxlist[0]},\$p" $FILE | grep '^PAIR' | head -n 1))
lenndx=${#content[*]}

if (( $lenndx == 0 )); then { echo "Fatal Error: invalid input file < $FILE >"; exit 0; } fi

# ndx: the index number of MAE value in echo PAIR line
for ((ndx=0; $ndx < $lenndx; ndx++))
do
    if [[ ${content[$ndx]} == $KEYWORD ]]; then break; fi
done
((ndx++))

if (($ndx >= $lenndx)); then { echo "Fatal Error: No KEYWORD=$KEYWORD is found"; exit 1; } fi


data=''
content=($(sed -n "${ndxlist[0]},\$p" $FILE | grep '^PAIR'))
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


# Process data: only MAE <= 0.2 will be left
# refstr  :  index number to be removed
# Tip: let < cnt > starting at zero
cnt=0
refstr=''
for i in $data
do
    ((cnt++))
    if [[ $i == 'nan' ]]; then { refstr="$refstr $cnt"; continue; } fi

    tmp=$(echo "0.2 - $i" | bc -l | grep '-')
    if [[ -n "$tmp" ]]; then refstr="$refstr $cnt"; fi
done


# Add Checking
tmp=$(sed -n "${ndxlist[0]}p" $FILE)
tmp=${tmp//*_}
tmp=${tmp%.*}
tmp=${tmp:=10000}
chklist=($tmp)

# Each Training's gennm
nmlist=()
for ((i=1; $i<${#ndxlist[*]}; i++))
do
    ((j=$i-1))
    vi=${ndxlist[$i]}
    vj=${ndxlist[$j]}

    nm=$(sed -n "$vj,${vi}p" $FILE | grep '^PAIR' | wc -l)
    nmlist+=($nm)

    tmp=$(sed -n "${vi}p" $FILE)
    tmp=${tmp//*_}
    tmp=${tmp%.*}
    tmp=${tmp:=10000}
    chklist+=($tmp)
done

# Checking < gennm >
# Because ndxlist contains the total number of lines of the FILE, remove it
((j = ${#chklist[*]} - 1))
for ((cnt=0; $cnt < $j; cnt++))
do
    ((cmp=$cnt+1))
    if [[ ${chklist[$cnt]} != $cmp ]]
    then
        echo "Fatal Error: < gennm > is not correctly identified in file < $FILE >"
        echo "             $(sed -n "${ndxlist[$cnt]}p" $FILE)"
        exit 1
    fi
done


# Now accumulating to Training number
tot=${nmlist[0]}
trylist=($tot)
for ((i=1; $i<${#nmlist[*]}; i++))
do
    ((tot=$tot+${nmlist[i]}))
    trylist+=($tot)
done


# Now based on refstr and trylist
# determine which Training and which number should be removed
# format of rmlist: "Training-nm PAIR-nm"
rmlist=()
for i in $refstr
do
    for ((j=0; $j<${#trylist[*]}; j++))
    do
        vj=${trylist[$j]}
        if (( $i <= $vj )); then break; fi
    done

    if (( $j == 0 ))
    then
        rmlist+=("1 $i")
    else
        ((tmp=$i-${trylist[$j-1]}))

        # Take care of number starting index
        ((ct=$j+1))
        rmlist+=("$ct $tmp")
    fi
done


# Now Do the cleaning
echo ''
for ((i=0; $i < ${#rmlist[*]}; i++))
do
    ndx=(${rmlist[$i]})
    x1=${ndx[0]}
    x2=${ndx[1]}

    if [[ -d Training_$x1 ]]
    then
        if (( $(ls Training_$x1/ | wc -l) <= 0 )); then continue; fi
        cd Training_$x1
        ftmp=Gas_$x2 
        if [[ -d $ftmp ]]; then { echo "Note: Deep Processing Training_$x1 Removing $ftmp"; rm -rf $ftmp; } fi

        ftmp=Liquid_$x2 
        if [[ -d $ftmp ]]; then { echo "Note: Deep Processing Training_$x1 Removing $ftmp"; rm -rf $ftmp; } fi

        ftmp=FEP_$x2 
        if [[ -d $ftmp ]]; then { echo "Note: Deep Processing Training_$x1 Removing $ftmp"; rm -rf $ftmp; } fi
        cd ../
    fi
done
echo ''


echo "Note: Cleaning current path"
rm -f *.err *.out
rm -f input_d.inp input_s.inpt


if [[ -d Test ]]
then
    echo "Note: Cleaning Test/ Folder"
    cd Test
    rm -rf Equil_*
    rm -f mdout* \#mdout*
    cd ../
fi


echo "Note: Cleaning Topfiles"
for fdir in $(ls | grep "Topfile_")
do
    if ! [[ -d $fdir ]]; then continue; fi
    cp -f $fdir/GAML_* .
    break
done
rm -rf Topfile_*
echo ''


for td in $(ls | grep 'Training_')
do
    if ! [[ -d $td ]]; then continue; fi
    if (( $(ls $td/ | wc -l) <= 0 ))
    then
        echo "Note: Removing Empty Folder $td"
        rmdir $td 
        continue
    fi
    cd $td
    for fdir in $(ls | grep "Gas_")
    do
        if ! [[ -d $fdir ]]; then continue; fi
        if (( $(ls $fdir/ | wc -l) <= 0 )); then continue; fi
        echo "Note: Processing $td Cleaning $fdir"
        cd $fdir
        rm -f mdout* \#mdout* energy*
        rm -f step* \#step*
        rm -f gas_min*
        rm -f gas_prod.cpt gas_prod.log gas_prod.trr gas_prod.edr
        cd ../
    done

    for fdir in $(ls | grep "Liquid_")
    do
        if ! [[ -d $fdir ]]; then continue; fi
        if (( $(ls $fdir/ | wc -l) <= 0 )); then continue; fi
        echo "Note: Processing $td Cleaning $fdir"
        cd $fdir
        rm -f mdout* \#mdout* energy*
        rm -f step* \#step*
        rm -f liq_min* liq_npt*
        rm -f liq_prod.cpt liq_prod.log liq_prod.trr
        cd ../
    done

    for fdir in $(ls | grep "FEP_")
    do
        if ! [[ -d $fdir ]]; then continue; fi
        if (( $(ls $fdir/ | wc -l) <= 0 )); then continue; fi
        echo "Note: Processing $fdir Cleaning:"
        cd $fdir
        rm -f bar_result.txt
        for dt in $(ls | grep "init_")
        do
            if ! [[ -d $dt ]]; then continue; fi
            echo -n "      $dt"
            cd $dt
            rm -f mdout* \#mdout* energy*
            rm -f step* \#step*
            rm -f grompp*
            rm -f fep_min* fep_npt*
            rm -f fep_prod.cpt fep_prod.log fep_prod.trr fep_prod.edr
            cd ../
        done
        echo ""
        cd ../
    done
    cd ../
    echo ''
done


if [[ -d originFile ]]
then
    echo "Note: Processing originFile"
    cd originFile
    rm -f grompp*
    rm -f *top *gro
    rm -rf Test/
    if [[ -d Training_origin ]]
    then
        cd Training_origin
        for i in $(ls)
        do
            if ! [[ -d $i ]]; then continue; fi
            cd $i
            rm -f mdout* \#mdout* step* \#step* energy*
            if [[ -n $(echo $i | grep "Gas") ]]
            then
                echo "Note: Processing originFile/Gas"
                rm -f gas_min*
                rm -f gas_prod.cpt gas_prod.edr gas_prod.log gas_prod.trr
            elif [[ -n $(echo $i | grep "Liquid") ]]
            then
                echo "Note: Processing originFile/Liquid"
                rm -f liq_min* liq_npt*
                rm -f liq_prod.cpt liq_prod.trr liq_prod.log
            elif [[ -n $(echo $i | grep "FEP") ]]
            then
                echo "Note: Processing originFile/FEP"
                rm -f bar_result.txt
                for j in $(ls)
                do
                    if ! [[ -d $j ]]; then continue; fi
                    cd $j
                    rm -f mdout* \#mdout* step* \#step* energy*
                    rm -f fep_min* fep_npt*
                    rm -f grompp*
                    rm -f fep_prod.log fep_prod.cpt fep_prod.trr
                    cd ../
                done
            else
                echo "Warning: Unknown Folder $i"
            fi
            cd ../
        done
        cd ../
    fi
    cd ../
fi

echo ''
echo "Everything is DONE"
