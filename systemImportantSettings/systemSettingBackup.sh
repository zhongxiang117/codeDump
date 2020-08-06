#!/bin/bash

# This script is used to backup system settings
# version 0.10 : July 31, 2020
# version 0.20 : more files can be added into < #MORE >

# update for more label
#MORE



VSCODE_PATH="/mnt/c/Users/11704/AppData/Roaming/Code/User"
VSCODE_PATH="/home/xiang/.config/Code/User"

VSCODE_PATH_FILE_SNIPPETS="${VSCODE_PATH}/snippets/globalSnippets.code-snippets"
VSCODE_PATH_FILE_KEYBINDINGS="${VSCODE_PATH}/keybindings.json"
VSCODE_PATH_FILE_SETTINGS="${VSCODE_PATH}/settings.json"
VSCODE_PATH_FILE_TASKS="${VSCODE_PATH}/tasks.json"

# folder for git backup
GIT_SYS_FOLDER="systemImportantSettings"

# exclusive for vscode
GIT_VSCODE_FOLDER="vsCodeSettingFiles"

GIT_URL=https://github.com/zhongxiang117/codeDump.git
GIT_FOLDER=${GIT_URL##*/}
GIT_FOLDER=${GIT_FOLDER%%.*}


CWD=$(pwd)
# this script
THIS_PATH_FILE="${CWD}/$0"



type git > /dev/null 2>&1 || { echo "Error: command < git > is not found"; exit 1; }
type uname > /dev/null 2>&1 || { echo "Error: command < uname > is not found"; exit 1; }

# decide machine type
MACHINE='linux'
if [[ -z $(uname -m | grep '86') ]]; then MACHINE='windows'; fi
echo "Note: identifying machine < $MACHINE >"


# check working path, git repository
if [[ ! -d ${GIT_FOLDER} ]]
then
    if [[ -z "$(git status > /dev/null 2>&1 | grep 'origin')" ]]
    then
        echo "Note: No repository is found, do you want to <git-clone> instead? y/yes, else quit"
        read tmp
        if [[ $tmp == 'y' || $tmp == 'yes' ]]
        then
            git clone ${GIT_URL}
        else
            echo "Note: you decided to quit..."
            echo "Note: exiting..."
            exit 1
        fi
    fi
fi


func_copyfiles() {
    # usage: copy $1 to current path if $1 not exist else copy to a "NEW-" name
    if (( $# != 1 )); then { echo "Error: check file path for <$1>"; exit 1; } fi
    fpath=$1
    if [[ -f ${fpath} ]]
    then
        filendx=${fpath##*/}
        if [[ -f ${filendx} ]]
        then
            if [[ -n "$(diff ${fpath} ${filendx})" ]]
            then
                newfilename="NEW-${filendx}"
                if [[ -f ${newfilename} ]]
                then
                    cnt=1
                    while true
                    do
                        if [[ -f "NEW-${cnt}-${filendx}" ]]; then ((cnt++)); else break; fi
                    done
                    newfilename="NEW-${cnt}-${filendx}"
                fi
                echo "Note: copying <$fpath> to <$newfilename>"
                cp ${fpath} ${newfilename}
            fi
        else
            echo "Note: copying <$fpath> to <$fpath>"
            cp ${fpath} .
        fi
    else
        echo "Warning: FILE < $fpath > not exist"
    fi
}



cd ${GIT_FOLDER}
if [[ ! -d ${GIT_SYS_FOLDER} ]]; then mkdir ${GIT_SYS_FOLDER}; fi
cd ${GIT_SYS_FOLDER}


# exclusive for this script
func_copyfiles ${THIS_PATH_FILE}

# exclusive for vscode
if [[ ! -d ${GIT_VSCODE_FOLDER} ]]; then mkdir ${GIT_VSCODE_FOLDER}; fi
cd ${GIT_VSCODE_FOLDER}
func_copyfiles ${VSCODE_PATH_FILE_KEYBINDINGS}
func_copyfiles ${VSCODE_PATH_FILE_SETTINGS}
func_copyfiles ${VSCODE_PATH_FILE_TASKS}
func_copyfiles ${VSCODE_PATH_FILE_SNIPPETS}
cd ../


# update for copying
#MORE



echo ''
echo ''
echo "Note: please checking < ${GIT_FOLDER}/${GIT_SYS_FOLDER} >"
echo "DONE everything"


