#!/bin/csh
alias HELP_argument 'eval \
echo "File < $0 > Help Usage:" \
echo "Command Line:" \
echo "    $0 [HEAD_File_Path] [Number_of_Simulation_Repeats]" \
echo "        Note: the sequence of input parameter is not important" \
echo "" \
echo "To see HEAD_File format, use:" \
echo "    $0 help HEAD" \
echo "Or" \
echo "    $0 HEAD help" \
echo ""\
echo ""\
echo "Input Parameter example:" \
echo "    1) with the HEAD existing, repeats numbr is 1" \
echo "        $0" \
echo "    2) with the HEAD existing, repeats numbr is rnm" \
echo "        $0 rnm" \
echo "    3) without the HEAD existing, repeats numbr is rnm" \
echo "        $0 rnm HEAD" \
echo "       or" \
echo "        $0 HEAD rnm" \
'

HELP_argument

echo done
