#!/bin/bash

# Files to be processed
FILE=MAE_PAIR_total.txt

# Number of entries to be shown
NUMBER=10

# Keyword to chosen
KEYWORD=%1


#
# END of user inputs
#

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Version 0.2   :   handle negative values
#   Version 0.3   :   smartly process "NaN"
#


FILE=${FILE:=MAE_PAIR_total.txt}
NUMBER=${NUMBER:=10}
KEYWORD=${KEYWORD:=\%1}


HELP_argument="Usage:    $0 [FILE] [NUMBER] [Keyword]

1)  FILE   : the input file needs to be processed. Default is < MAE_PAIR_total.txt >
2)  NUMBER : the printout number of entries. Default is < 10 >
3)  Keyword: only Four Keywords can be chosen: < %1 %2 %3 %4 >
            < %1 >: MAE (default);  < %2 >: Density; < %3 >: HVAP;  < %4 >: FEP
"

nmargv=$#
if (( $nmargv == 0 || $nmargv == 2 || $nmargv == 3))
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



if (( $nmargv == 0 ))
then
    if ! [[ -f $FILE ]]
    then
        echo "Error: < $FILE > file does not exist, quitting..."
        echo ''
        exit 1
    fi
elif (( $nmargv == 1 ))
then
    if [[ -f $1 ]]
    then
        FILE=$1
    elif [[ $1 =~ ^[0-9]+$ ]]
    then
        NUMBER=$1
    elif [[ $1 =~ %[1,2,3,4]+ ]]
    then
        KEYWORD=$1
    else
        echo "Error: invalid input parameter < $1 >"
        exit 1
    fi
elif (( $nmargv == 2 ))
then
    if [[ -f $1 ]]
    then
        FILE=$1
        if [[ $2 =~ ^[0-9]+$ ]]
        then
            NUMBER=$2
        elif [[ $2 =~ %[1,2,3,4]+ ]]
        then
            KEYWORD=$2
        else
            echo "Error: invalid input parameter < $2 >"
        fi
    elif [[ -f $2 ]]
    then
        FILE=$2
        if [[ $1 =~ ^[0-9]+$ ]]
        then
            NUMBER=$1
        elif [[ $1 =~ %[1,2,3,4]+ ]]
        then
            KEYWORD=$1
        else
            echo "Error: invalid input parameter < $1 >"
        fi
    elif [[ -f $FILE ]]
    then
        if [[ $1 =~ ^[0-9]+$ && $2 =~ %[1,2,3,4]+ ]]
        then
            NUMBER=$1
            KEYWORD=$2
        elif [[ $2 =~ ^[0-9]+$ && $1 =~ %[1,2,3,4]+ ]]
        then
            NUMBER=$2
            KEYWORD=$1
        else
            echo "Error: invalid input parameter < $1 > or < $2 >"
            exit 1
        fi
    else
        echo "Error: invalid input parameter < $1 > or < $2 >"
        exit 1
    fi
else
    if [[ -f $1 ]]
    then
        FILE=$1
        if [[ $2 =~ ^[0-9]+$ && $3 =~ %[1,2,3,4]+ ]]
        then
            NUMBER=$2
            KEYWORD=$3
        elif [[ $3 =~ ^[0-9]+$ && $2 =~ %[1,2,3,4]+ ]]
        then
            NUMBER=$3
            KEYWORD=$2
        else
            echo "Error: invalid input parameter < $2 > or < $3 >"
            exit 1
        fi
    elif [[ -f $2 ]]
    then
        FILE=$2
        if [[ $1 =~ ^[0-9]+$ && $3 =~ %[1,2,3,4]+ ]]
        then
            NUMBER=$1
            KEYWORD=$3
        elif [[ $3 =~ ^[0-9]+$ && $1 =~ %[1,2,3,4]+ ]]
        then
            NUMBER=$3
            KEYWORD=$1
        else
            echo "Error: invalid input parameter < $1 > or < $3 >"
            exit 1
        fi
    elif [[ -f $3 ]]
    then
        FILE=$3
        if [[ $1 =~ ^[0-9]+$ && $2 =~ %[1,2,3,4]+ ]]
        then
            NUMBER=$1
            KEYWORD=$2
        elif [[ $2 =~ ^[0-9]+$ && $1 =~ %[1,2,3,4]+ ]]
        then
            NUMBER=$2
            KEYWORD=$1
        else
            echo "Error: invalid input parameter < $1 > or < $2 >"
            exit 1
        fi
    else
        echo "Error: invalid input parameter < $1 > or < $2 > or < $3 >"
        exit 1
    fi
fi


# transfter KEYWORD to real KEYWORD
if [[ $KEYWORD == %1 ]]; then KEYWORD=MAE; elif [[ $KEYWORD == %2 ]]; then KEYWORD=Density; \
elif [[ $KEYWORD == %3 ]]; then KEYWORD=HVAP; else KEYWORD=FEP; fi


echo "Note: Setting processing file to < $FILE >..."
echo "Note: Setting number of printouts to < $NUMBER >..."
echo "Note: Setting Keyword to < $KEYWORD >..."
echo ""


content=($(cat $FILE | grep '^PAIR' | head -n 1))
lenndx=${#content[*]}

if (( $lenndx == 0 )); then { echo "Error: invalid input file < $FILE >"; exit 0; } fi

for ((ndx=0; $ndx < $lenndx; ndx++))
do
    if [[ ${content[$ndx]} == $KEYWORD ]]; then break; fi
done
((ndx++))


content=($(cat $FILE | grep '^PAIR'))
lentot=${#content[*]}

if (( $lentot == 0 )); then { echo "Error: no inputs after processing for < $FILE >"; exit 0; } fi

data=''
for ((i=$ndx; $i<$lentot; i=$i+$lenndx))
do
    data="$data ${content[$i]}"
done

# Annoying thing is `sort` cannot handle negative numbers
# We will use an different way
#dataFilter=$(echo $data | tr ' ' '\n' | sort -n | head -n $NUMBER)

# Convert data to an arrya
dataArray=($data)


# Bubble sort -- Very Very Very time consuming 
for ((i=1; $i < ${#dataArray[*]}; i++))
do
    for ((j=0; $j < $i; j++))
    do
        vi=$(echo "${dataArray[$i]}" | sed 's/-//' | tr '[A-Z]' '[a-z]')
        vj=$(echo "${dataArray[$j]}" | sed 's/-//' | tr '[A-Z]' '[a-z]')

        if [[ $vi == 'nan' || $vj == 'nan' ]]
        then
            toi='nan'
            if [[ $vj == 'nan' ]]
            then
                toj=${dataArray[$i]}
            else
                toj=${dataArray[$j]}
            fi
            dataArray[$i]=$toi
            dataArray[$j]=$toj
        else
            err=$(echo "$vi - $vj" | bc -l | grep '-')
            if [[ -n "$err" ]]
            then
                tmp=${dataArray[$i]}
                dataArray[$i]=${dataArray[$j]}
                dataArray[$j]=$tmp
            fi
        fi
    done
done


dataFilter=''
for ((i=0; $i<$NUMBER; i++))
do
    v=${dataArray[$i]}
    dataFilter="$dataFilter $v"
    echo "No. $cnt:   $KEYWORD  $v"
done
echo ""


echo "Do you want to see full entries? y/yes, else not"
read t
if [[ -n "$t" && ( $t == 'y' || $t == 'yes' ) ]]
then
    # select all corrected index but will result in an random sequence
    maestr=''
    refstr=''
    cnt=0
    for i in $data
    do
        for j in $dataFilter
        do
            if [[ $i == $j ]]
            then
                refstr="$refstr $cnt"
                maestr="$maestr $i"
                break
            fi
        done
        ((cnt++))
    done

    pairlist=()
    for ((i=$lenndx; $i<=$lentot; i=$i+$lenndx))
    do
        line=''
        for ((j=$i-$lenndx; $j<$i; j++))
        do
            line="$line ${content[$j]}"
        done
        pairlist+=("$line")
    done


    filterlist=()
    for i in $refstr
    do
        filterlist+=("${pairlist[$i]}")
    done

    # convert refstr to reflist
    reflist=($refstr)


    # avoid any same mae value repeats
    # Idea is using an array holding USED index of $dataFilter
    # everytime encouter identify, check this array first
    cnt=0
    replist=()
    for i in $dataFilter
    do
        t=0
        for j in $maestr
        do
            if [[ $i == $j ]]
            then
                bo=true
                for ((k=0; $k<${#replist[*]}; k++))
                do
                    if (( $cnt == ${replist[$k]} )); then { bo=false; break; } fi
                done

                # update replist
                replist+=($cnt)

                if $bo
                then
                    tmp=${reflist[$t]}
                    echo "${pairlist[$tmp]}" | sed 's/^ //' | tr ' ' '\t'
                    break
                fi
            fi
            ((t++))
        done
        ((cnt++))
    done
fi
