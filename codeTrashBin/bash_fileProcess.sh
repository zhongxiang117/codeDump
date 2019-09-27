#!/bin/bash

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   version 0.2: process files, mkdir; and rename them to MOL.*
#
#

FILE=files



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

content=''
for f in $(cat $FILE)
do
    if [[ -f $f ]]
    then
        name=${f%.*}
        content="$content $name"
        ext=${f/*.}

        if [[ -z "$ext" ]]; then { echo "Error: wrong file extension"; continue; } fi

        echo "Note: Processing $f"
        sed -i 's/UNK/MOL/' $f

        if [[ $ext == pdb ]]
        then
            sed -i '/^CONECT/d; /^TER/d' $f
        fi

        if [[ $ext == itp ]]
        then
            nm=$(sed -n '/\[ pairs \]/=' $f)
            if [[ -n "$nm" ]]
            then
                sed -i "$nm,\$d" $f
            fi

            # double check
            while read -r line
            do
                pl=($line)
                if (( ${#pl[*]} == 3 ))
                then
                    bo=true
                    for i in ${pl[*]}
                    do
                        if ! [[ $i =~ ^[0-9]+$ ]]
                        then
                            bo=false
                            break
                        fi
                    done
                    if $bo
                    then
                        echo "Error: on processing $f"
                        echo ''
                        echo "Error Line: $line"
                        echo ''
                        break
                    fi
                fi
            done < $f
        fi
    else
        echo "Warning: $f does not exist"
    fi
done

if [[ -z "$content" ]]; then { echo "Error: no inputs"; exit 1; } fi

content=($content)
unicontent="${content[0]}"
for ((i=1; $i<${#content[*]}; i++))
do
    bo=true
    for ((j=0; $j<$i; j++))
    do
        if [[ ${content[$j]} == ${content[i]} ]]
        then
            bo=false
            break
        fi
    done
    if $bo; then unicontent="$unicontent ${content[$i]}"; fi
done

# move files to its directory
for fnm in $unicontent
do
    if ! [[ -d $fnm ]]
    then
        mkdir $fnm
        mv $fnm\.* $fnm/
    fi
done

# rename file name
for dir in $unicontent
do
    cd $dir
    if [[ -f $dir.itp ]]
    then
        if ! [[ -f MOL.itp ]]; then mv $dir.itp MOL.itp; fi
    fi

    if [[ -f $dir.pdb ]]
    then
        if ! [[ -f MOL.pdb ]]; then mv $dir.pdb MOL.pdb; fi
    fi

    if [[ -f $dir.gro ]]
    then
        if ! [[ -f MOL.gro ]]; then mv $dir.gro MOL.gro; fi
    fi
    cd ../

	echo "$dir" >> NEWfiles
done




echo ''
echo "Done Everything"
