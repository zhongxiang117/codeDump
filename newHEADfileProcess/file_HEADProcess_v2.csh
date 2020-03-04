#!/bin/csh

set nmargv =  ${#argv}

alias HELP_argument 'eval \
echo "CSHFILE < $0 > Help Usage:" \
echo "Command Line:" \
echo "    $0 [HEAD_File_Path] [Number_of_Simulation_Repeats]" \
echo "        <Note: the sequence of input parameter is not important>" \
echo "" \
echo "To see HEAD_File format, use:" \
echo "    $0 help HEAD" \
echo "Or" \
echo "    $0 HEAD help" \
echo "" \
echo "" \
echo "Input Parameter example:" \
echo "    1) with HEAD file existing, repeats numbr is 1" \
echo "        $0" \
echo "    2) with HEAD file existing, repeats numbr is rnm" \
echo "        $0 rnm" \
echo "    3) without HEAD file existing, repeats numbr is rnm" \
echo "        $0 rnm HEAD" \
echo "       or" \
echo "        $0 HEAD rnm" \
'

alias HELP_format  'eval \
echo "CSHFILE < $0 > Help Usage:" \
echo "Command Line:" \
echo "    $0 [HEAD_File_Path] [Number_of_Simulation_Repeats]" \
echo "        <Note: the sequence of input parameter is not important>" \
echo "" \
echo "To see Command line input format, use:" \
echo "    $0 help       OR:   $0 -h" \
echo "OR: $0 --help     OR:   $0 -help" \
echo "" \
echo "" \
echo "Input HEAD file format:" \
echo "" \
echo "The first six lines have to be in sequence with 12 parameters" \
echo "" \
echo "Line-1:    G09          [any-string]" \
echo "Line-2:    GTYPE        [any-string]" \
echo "Line-3:    TOTCON       [integer_numbr]" \
echo "Line-4:    GUPCMD       [integer_numbr]" \
echo "Line-5:    QM_COUNT     [integer_numbr]" \
echo "Line-6:    FREQUENCY    [integer_numbr]" \
'

if ( $nmargv == 0 ) then
    :
else if ( $nmargv == 1 ) then
    if ( "$1" == 'help' || "$1" == '-h' || "$1" == '--help' || \
         "$1" == '-help' ) then
        HELP_argument
        exit 0
    endif
else if ( $nmargv == 2 ) then
    if ( ( "$1" == 'help' && "$2" == 'HEAD' ) || \
         ( "$1" == 'HEAD' && "$2" == 'help' ) ) then
        HELP_format
        exit 0
    endif
else
    echo 'Error: Number of input parameters were wrong, quitting...'
    echo ''
    exit 1
endif



if ( $nmargv == 0 ) then
    echo 'Note: Setting number of simulation repeats to < 1 >...'
    set repnm = 1
    if ! ( -f HEAD ) then
        echo 'Error: < HEAD > file does not exist, quitting...'
        echo ''
        exit 1
    endif
    echo 'Note: Using < HEAD > file as processing file...'
    set headfile = 'HEAD'
else if ( $nmargv == 1 ) then
    if ( -f ${argv[1]} ) then
        echo 'Note: Setting number of simulation repeats to < 1 >...'
        set repnm = 1
        echo "Note: Setting < ${argv[1]}> as processing file..."
        set headfile = ${argv[1]}
    else
        set tmp = `echo ${argv[1]} | grep '^[0-9]\+$'`
        if ( $tmp == '' ) then
            echo "Error: The input parameter < ${argv[1]} > was wrong, quitting..."
            echo ''
            exit 1
        endif
        echo "Note: Setting number of simulation repeats to < $tmp >..."
        set repnm = $tmp
        if ! ( -f HEAD ) then
            echo 'Error: < HEAD > file does not exist, quitting...'
            echo ''
            exit 1
        endif
        echo 'Note: Using < HEAD > file as processing file...'
        set headfile = 'HEAD'
    endif
else
    if ( -f ${argv[1]} ) then
        set tmp = `echo ${argv[2]} | grep '^[0-9]\+$'`
        if ( $tmp == '' ) then
            echo "Error: The input parameter < ${argv[2]} > was wrong, quitting..."
            echo ''
            exit 1
        endif
        echo "Note: Setting number of simulation repeats to < $tmp >..."
        set repnm = $tmp
        echo "Note: Setting < ${argv[1]} > as processing file..."
        set headfile = "${argv[1]}"
    else
        set tmp = `echo ${argv[1]} | grep '^[0-9]\+$'`
        if ( $tmp == '' ) then
            echo "Error: The input parameter < ${argv[1]} > was wrong, quitting..."
            echo ''
            exit 1
        endif
        echo "Note: Setting number of simulation repeats to < $tmp >..."
        set repnm = $tmp
        if ( -f ${argv[2]} ) then
            echo "Note: Setting < ${argv[2]} > as processing file..."
            set headfile = "${argv[2]}"
        else
            echo "Error: The input parameter < ${argv[2]} > was wrong, quitting..."
            echo ''
            exit 1
        endif
    endif
endif


set contents = `head -n 6 $headfile`
if ( ${#contents} != 12 || ${contents[1]} != 'G09' || \
     ${contents[3]} != 'GTYPE' || ${contents[5]} != 'TOTCON' || \
     ${contents[7]} != 'GUPCMD' || ${contents[9]} != 'QM_COUNT' || \
     ${contents[11]} != 'FREQUENCY' ) then
    echo "Error: The format for input < $headfile > file is wrong!"
    echo ''
    exit 1
endif


set totcon = `echo ${contents[6]} | grep '^[0-9]\+$'`
set fgup = `echo ${contents[8]} | grep '^[0-9]\+$'`
set qmcnt = `echo ${contents[10]} | grep '^[0-9]\+$'`
set fcnt = `echo ${contents[12]} | grep '^[0-9]\+$'`
if ( "$totcon" == '' || "$fgup" == '' || "$qmcnt" == '' || "$fcnt" == '' || \
      $totcon <= 0 || $fgup <= 0 || $qmcnt <= 1 || $fcnt < 1 ) then
    echo "Error: The format for input < $headfile > file is wrong!"
    echo ''
    exit 1
endif


if ( $fcnt - 1 == 0 ) then
    echo 'Error: Currently only the type of UpdateGaussian is a Number is supported'
    echo 'Error: For more information, please refer IFLINE file'
    echo ''
    exit 1
endif


set bool_broken = 'n'
if ( `expr \( $qmcnt - 1 \) % 3` != 0 ) set bool_broken = 'y'
@ rcnt = ( $qmcnt - 1 ) / 3
set Ygauss = `grep GAUSSACCEPT $headfile | wc -l`
set Yaenet = `grep AENETACCEPT $headfile | wc -l`
set Yaeg09 = `grep AEG09ACCEPT $headfile | wc -l`
set Yginit = `grep 'INIT DONE' $headfile | wc -l`
if ( $Ygauss + $Yaenet > $rcnt || $Yaeg09 > $Yaenet || $totcon < $rcnt - 1 ) then
    echo "Error: The format for input < $headifile > file is wrong!"
    echo ''
    exit 1
endif

if ( `expr $Yginit % 3` != 0 ) set bool_broken = 'y'
if ( `expr $rcnt % $repnm` != 0 ) set bool_broken = 'y'
if ( `expr $totcon % \( $fcnt - 1 \)` != 0 ) set bool_broken = 'y'
if ( `expr $rcnt % $repnm` != 0 ) set bool_broken = 'y'



@ totgauss = $fcnt + $Yaeg09 + $Yginit
@ totaenet = $rcnt - $fcnt

set ratio_node = `echo "$Ygauss / $totgauss * 100" | bc -l`
set ratio_node = `printf '%.2f' $ratio_node`
set ratio_node = "$ratio_node""%"

set ratio_aenet = `echo "$Yaenet / $totaenet * 100" | bc -l`
set ratio_aenet = `printf '%.2f' $ratio_aenet`
set ratio_aenet = "$ratio_aenet""%"

set ratio_aeg09 = `echo "$Yaeg09 / $Yaenet * 100" | bc -l`
set ratio_aeg09 = `printf '%.2f' $ratio_aeg09`
set ratio_aeg09 = "$ratio_aeg09""%"

set ratio_gaund = `echo "$Ygauss / $fcnt * 100" | bc -l`
set ratio_gaund = `printf '%.2f' $ratio_gaund`
set ratio_gaund = "$ratio_gaund""%"

set ratio_gauss = `echo "$Yaeg09 / $totgauss * 100" | bc -l`
set ratio_gauss = `printf '%.2f' $ratio_gauss`
set ratio_gauss = "$ratio_gauss""%"

echo ''
echo "For input simulation:"
if ( $bool_broken == 'y' ) echo 'Warning: the simulation is broken or not yet finished'
echo ''
echo "Total configuration          = $totcon"
if ( $bool_broken == 'n' ) then
    @ xf = $totcon / ( $rcnt - 1 ) / $repnm
    echo "Solute Move Frequency        = $xf"
endif
echo "Update Gaussian Frequency    = $fgup"
echo "Number of simulation Repeats = $repnm"
echo ''
echo "Total Number of Gaussian Calc         = $totgauss"
echo "Total Number of Init Gaussian Calc    = $Yginit"
echo "Total Number of Node Gaussian Calc    = $fcnt"
echo "Total Number of Node AENET Calc       = $totaenet"
echo "Total Number of Tossing Gaussian Calc = $Yaeg09"
echo ''
echo "The Accept Ratio Between Node Gaussian and Total Gaussian    = $ratio_node"
echo "The Accept Ratio for Node Gaussian in routine Node Calc      = $ratio_gaund"
echo "The Accept Ratio Between Node AENET and Total AENET          = $ratio_aenet"
echo "The Accept Ratio Between Tossing Gaussian and Accept AENET   = $ratio_aeg09"
echo "The Accept Ratio Between Tossing Gaussian and Total Gaussian = $ratio_gauss"
echo ''

