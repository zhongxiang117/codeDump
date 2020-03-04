#!/bin/csh

echo $0
echo good
exit 0

if ( $#argv == 1 ) then
    if ( "$1" == '-h' || "$1" == '-help' || \
         "$1" == '--help' || "$1" == 'help' ) then
            echo good
    endif
else if ( $#argv == 2 ) then
    switch ($argv[1])
        case '-h':
        case '-help':
        case '--help':
        case 'help':
            set par = "$2"
            breaksw
        
        case 'HEAD':
        case 'File':
            echo good
            breaksw

          default:
              breaksw
      endsw
endif





echo done





