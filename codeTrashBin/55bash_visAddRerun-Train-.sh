#!/bin/bash
# -*- coding: utf-8 -*-

#BSUB -n 1
#BSUB -q gpu
#BSUB -W 168:00
#BSUB -J IL
#BSUB -P prmt
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -x

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# version 0.10
# version 0.20  More robust
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# global parameter
bool_check=true


# processing
bool_process=false
# unit in ps
analysis_beginning_time=20


# Folder needs rerun
# only the largest number can be auto-identified
TOPFOLDER=Topfile_13-1

# which subfolder to add, integer array, string
DIRVISNM="57 58 60 61 63 65 66 67 68 69"

# cos-acceleration, string
VISCOSNM="0.15 0.25 0.35 0.4"


gmx=gmx



#
# End of user inputs
#

gmx=${gmx:=gmx}
if $bool_check; then bool_process=false; fi
if $bool_process; then bool_check=false; fi


# auto-identify TOPFOLDER
if [[ -z "$TOPFOLDER" ]]
then
    nm=""
    for f in $(ls | grep Topfile)
    do
        if [[ -d $f ]]
        then
            f=${f##*_}
            if [[ "$f" =~ ^[0-9]+$ ]]; then nm="$nm $f"; fi
        fi
    done
    nm=($(echo "$nm" | tr ' ' '\n' | sort -r -n))
    TOPFOLDER=Topfile_${nm[0]}
fi


GROMPP=grompp_prod_viscosity.mdp
if [[ ! $bool_process && ! -f $GROMPP ]]; then { echo "Error: < vis-grompp > does not exist"; exit 1; } fi

if [[ ! $bool_process && ! -d $TOPFOLDER ]]; then { echo "Error: < $TOPFOLDER > does not exist"; exit 1; } fi

MINI=MINI_$TOPFOLDER
if [[ ! $bool_process && ! -d $MINI ]]; then { echo "Error: Folder < $MINI > does not exist"; exit 1; } fi

TRAIN=Training_$TOPFOLDER
if [[ ! -d $TRAIN ]]; then { echo "Error: Folder < $TRAIN > does not exist"; exit 1; } fi



# double check DIRVISNM
if ! $bool_process
then
    # get training topnm
    topnm=""
    for i in $(ls $TOPFOLDER)
    do
        t=${i##*_}
        t=${t%%\.*}
        if [[ "$t" =~ ^[0-9]+$ ]]; then topnm="$topnm $t"; fi
    done
    topnm=$(echo "$topnm" | tr ' ' '\n' | sort -n)

    for i in $DIRVISNM
    do
            if [[ -z "$(echo $topnm | grep $i)" ]]
            then
                echo "Error: rerun number < $i > is not in < $TOPFOLDER >"
                exit 1
            fi

            if [[ ! ( -f $MINI/min_${i}.gro && -f $MINI/min_${i}.log ) ]]
            then
                echo "Error: rerun number < $i > does not have < mini-gro or mini-log > file"
                exit 1
            fi

            if [[ ! -f $TRAIN/dir_$i/prod.gro ]]
            then
                echo "Error: rerun number < $i > does not have < prod-gro > file"
                exit 1
            fi

            for vnm in $VISCOSNM
            do
                if [[ -d $TRAIN/dir_$i/vis_$vnm ]]
                then
                    echo "Error: rerun number < $i > < vis_$vnm > already exisit"
                    exit 1
                fi
            done
    done
fi


# now get forces and exit for bool_check
if $bool_check
then
    for i in $DIRVISNM
    do
        t=$(grep 'Norm of force' $MINI/min_${i}.log)
        echo "dir_$i     $t"
    done
    echo ""
    echo "Will process < $TOPFOLDER >"
    echo "Folder modify < $DIRVISNM >"
    echo "Will add vis-nm < $VISCOSNM >"
    echo "Done checking"

    exit 0
fi


if $bool_process
then
    cd $TRAIN
    pwd

    echo '' > rstfile_vis
    for fd in $DIRVISNM
    do
        if [[ ! -d dir_$fd ]]; then continue; fi
        cd dir_$fd
        for vf in $VISCOSNM
        do
            if [[ -d vis_$vf ]]
            then
                cd vis_$vf
                if [[ -f vis.gro && -f vis.edr ]]
                then
                    echo "@@VIS ${fd} ${vf}" >> ../../rstfile_vis
                    echo {36..42} | $gmx energy -f vis.edr -b $analysis_beginning_time >> ../../rstfile_vis
                    echo "" >> ../../rstfile_vis

                    # cleaning
                    rm -f energy* \#energy* mdout* \#mdout*
                    rm -f step* \#step*
                    rm -f vis.log vis.trr vis_prev.cpt vis.cpt
                    rm -f core*
                else
                    echo "${vf} was broken" >> ../../rstfile_vis
                    echo '' >> ../../rstfile_vis
                    rm -f *
                fi

                cd ../
            else
                echo "${vf} not exist" >> ../rstfile_vis
                echo '' >> ../rstfile_vis
                echo '' >> ../rstfile_vis
            fi
        done
        cd ../
    done


    echo "Note: Processing..."
    echo ''
    echo ''

    # format: "dirnm   visnm   1/visValue": "57   0.15   10.3452"
    dvislist=()

    ## nmlist: hold line numbers for all "@@VIS"
    nmlist=()
    for i in $(grep -n '@@VIS' rstfile_vis)
    do
        if [[ -n "$(echo $i | grep ':')" ]]; then nmlist+=(${i%%:*}); fi
    done
    nmlist+=($(cat rstfile_vis | wc -l))

    for ((i=1; $i<${#nmlist[*]}; i++))
    do
        ((j=$i-1))
        beg=${nmlist[j]}
        end=${nmlist[i]}
        visrep=($(sed -n "$beg,${end}p" rstfile_vis | grep 1/Visco))
        visrep=${visrep[1]}
        if [[ -z "$visrep" ]]; then visrep=NaN; fi

        nm=($(sed -n "${beg}p" rstfile_vis))
        dirnm=${nm[1]}
        visnm=${nm[2]}

        dvislist+=("$dirnm $visnm $visrep")
    done

    # using "dirnm" collect all other values
    i=0
    while ((i<${#dvislist[*]}))
    do
        ilist=(${dvislist[$i]})
        idir=${ilist[0]}

        vcol="vis_${ilist[1]}  ${ilist[2]}"

        for ((j=$i+1; $j<${#dvislist[*]}; j++))
        do
            jlist=(${dvislist[$j]})
            jdir=${jlist[0]}
            
            if [[ $jdir == $idir ]]
            then
                jvalue="vis_${jlist[1]}  ${jlist[2]}"
                vcol="$vcol    $jvalue"
            else
                break
            fi
        done

        i=$j

        echo "$idir   $vcol"
    done

    echo ''
    echo ''
    echo "Done-processing for < $TOPFOLDER >"
    exit 0
fi


# Training
cd $TRAIN
for fd in $DIRVISNM
do
    cd dir_$fd

    for vnm in $VISCOSNM
    do
        mkdir vis_$vnm
        cd vis_$vnm

        mdpfile=grompp_prod_vis_${vnm}.mdp
        topfile=top_$fd.top

        cp ../../../$GROMPP $mdpfile

        sed -i "s/VAR_COS/$vnm/" $mdpfile

        grompp -f $mdpfile -c ../prod.gro -p ../../../$TOPFOLDER/$topfile -o vis.tpr

        if [[ -f vis.tpr ]]
        then
            $gmx mdrun -deffnm vis
            wait
        else
            echo "Error: on < dir_$fd > < vis_$vnm.tpr >"
        fi

        cd ../
    done

    cd ../
done







