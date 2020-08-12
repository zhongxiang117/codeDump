# backup for import information

# windows-ubuntu root path
alias ubuntu-root-abspath="echo 'C:\Users\USER\AppData\Local\Packages\CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc\LocalState\rootfs'"


# Byte Order Mark for Python
# PEP 263 -- Defining Python Source Code Encodings
alias byte-order-mark="echo '# -*- coding: utf-8 -*-'"


# exclusively for Gaussian 16
#g16root="/home/USER/Apps/Gaussian16"
#    GAUSS16_SCRDIR="/home/USER/Temp/ScratchG16"
#    export g16root GAUSS16_SCRDIR
#    . $g16root/g16/bsd/g16.profile


# For wsl, set Color like ubuntu theme
#
#    Slot 1:  Red: 78,  Green: 154, Blue: 6  -- change it to < #3 > to avoid < vim > color showing problem
#    Slot 2:  Red: 52,  Green: 101, Blue: 164
#    Slot 3:  Red: 48,  Green: 10,  Blue: 36
#    Slot 4:  Red: 6,   Green: 152, Blue: 154
#    Slot 5:  Red: 204, Green: 0,   Blue: 0
#    Slot 6:  Red: 117, Green: 80,  Blue: 123
#    Slot 7:  Red: 196, Green: 160, Blue: 0
#    Slot 8:  Red: 211, Green: 215, Blue: 207
#    Slot 9:  Red: 85,  Green: 87,  Blue: 83
#    Slot 10: Red: 114, Green: 159, Blue: 207
#    Slot 11: Red: 138, Green: 226, Blue: 52
#    Slot 12: Red: 52,  Green: 226, Blue: 226
#    Slot 13: Red: 239, Green: 41,  Blue: 41
#    Slot 14: Red: 173, Green: 127, Blue: 168
#    Slot 15: Red: 252, Green: 233, Blue: 79
#    Slot 16: Red: 238, Green: 238, Blue: 238
#
# Then for color tab
#
#    Set Text and Popup Text to < #16 >
#    Set Background and Popup Background to < #3 >
#
# Finally, install and set ubuntu Mono font




# REG_SZ          : String
# REG_MULTI_SZ    : Multi-String
# REG_EXPAND_SZ   : Expandable Variables <%>
# REG_BINARY      : Binary
# REG_DWORD       : <32-bit> -- decimal & hexdecimal
# REG_QWORD       : <64-bit> -- decimal & hexdecimal
#
# HKEY_CLASSES_ROOT     HKCR    Stores file association and COM object registration
# HKEY_CURRENT_USER     HKCU    Stores data associated with the account currently logged on
# HKEY_LOCAL_MACHINE    HKLM    Stores system-related information
# HKEY_USERS            HKU     Stores  information about all the accounts on the machine
# HKEY_CURRENT_CONFIG   HKCC    Stores information about the current machine profile
#
# Open With context menu.
#
#   Computer\HKEY_CLASSES_ROOT\*\shell\Open with wordpad\Command\KEY_PATH_TO_EXE
#   Computer\HKEY_CURRENT_USER\Software\Microsoft\Windows\CurrentVersion\Explorer\FileExts\.sh
#
# Desktop Context Menu
#   HKEY_CLASSES_ROOT\Directory\Background\shell
#
# < OpenWithList > Variables under:
#   HKEY_CLASSES_ROOT/Applications
#
# For all of them, the extend shell (via <K_Ctrl>) is with < shellex > instead
#
#
# HKEY_CLASSES_ROOT\*\OpenWithList
# HKEY_CLASSES_ROOT\*\shellex\ContextMenuHandlers
# HKEY_CLASSES_ROOT\Directory\Background
# HKEY_CLASSES_ROOT\Directory\shell
# HKEY_CLASSES_ROOT\Directory\shellex\ContextMenuHandlers
# HKEY_CLASSES_ROOT\Drive\shell
# HKEY_CLASSES_ROOT\Drive\shellex\ContextMenuHandlers
# HKEY_CLASSES_ROOT\Folder\shell
# HKEY_CLASSES_ROOT\Folder\shellex\ContextMenuHandlers 
#
# Tree:
#
# HKCR/*/shell/
#   CustomOpen/
#       Command/    --> REG_SZ{KEY_PATH_TO_EXE,"%1"}
#       Icon        --> REG_SZ{.exe,.ico}
#       Position    --> REG_SZ{top,bottom}


# Remove SkyDrive Pro in Context menu
#   Under      < HKEY_CLASSES_ROOT\AllFilesystemObjects\shell >
#   Delete key < SPFS.ContextMenu >


# add git ssh-key
#
# cmd: ssh-keygen -t rsa -C "email@addr.com"
# default under < .ssh/ > will have files < id_rsa > (private) & < id_rsa.pub > (public)
# paste content in file < id_rsa.pub > to git account
#
# cmd: ssh-agent
#      $? --> 0 running with passpharse, 1 running without passpharse, 2 not running
# cmd: ssh-add


# restart ubuntu GUI
# ALT + F2  ->  r
# gnome-shell --replace & disown

# desktop file location
# /usr/share/applications
# ~/.local/share/applications




