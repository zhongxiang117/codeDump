#!/bin/bash


t=$(type code)
if [[ -z $t ]]; then echo "Error: VSCode is not found"; exit; fi


MY_EXTENSIONS_ID='
alefragnani.bookmarks
coenraads.bracket-pair-colorizer
auchenberg.vscode-browser-preview
ms-vscode.cpptools
ms-vscode.cmake-tools
formulahendry.code-runner
streetsidesoftware.code-spell-checker
eff-hykin.code-eol
msjsdiag.debugger-for-chrome
hediet.vscode-drawio
janisdd.vscode-edit-csv
grapecity.gc-excelviewer
letrieu.expand-region
mkxml.vscode-filesize
hansec.fortran-ls
sirtori.indenticator
lonefy.vscode-js-css-html-formatter
james-yu.latex-workshop
sissel.shopify-liquid
ritwickdey.liveserver
yzhang.markdown-all-in-one
yzane.markdown-pdf
gimly81.matlab
krvajalm.linter-gfortran
ms-python.python
lextudio.restructuredtext
killalau.vscode-liquid-snippets
rashwell.tcl
visualstudioexptteam.vscodeintellicode
ms-vscode.vs-keybindings
vscode-icons-team.vscode-icons
omoki1207.pdf
ms-toolsai.jupyter
'

# pip3 install jupyter notebook ipykernel


for eid in $MY_EXTENSIONS_ID; do
    echo "Note: code installing $eid ..."
    code --install-extension $eid
done



