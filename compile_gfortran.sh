###############################################################
#                                                             #
# Brief description:                                          #
# This script compile and link PSIDE source codes with        #
# gfortran compiler.                                          #
# The use cases is:                                           #
#   ./compile_src_gfortran.sh test_problem               #
# where test_problem is the fortran source file which descri- #
# bes the problem to be solved by PSIDE.                      # 
#                                                             #
# Author: Diego T. Volpatto                                   #
# Contact: volpatto@lncc.br                                   #
#                                                             #
###############################################################

# Check if problem is supplied
if [ -z "$1" ]; then
    echo "Error: The problem was not supplied"
    exit 1
fi

# Verifing gfortran version
echo "Compiler version check:"
echo -e "+++ Version:\n"
comando="gfortran --version"; eval $comando
if [ $? -ne 0 ]; then
    echo "+++ Compiler can not be checked"
    echo "+++ Status: Failed"
    exit 1
fi
echo -e "+++ Status: OK\n"

# Clearing
echo "Clearing binaries:"
binDir="bin/"
if [ -d "$binDir" ]; then
  comando="cd $binDir"; eval $comando
  comando="rm *.o"; eval $comando
else
  comando="mkdir bin"; eval $comando
  comando="cd $binDir"; eval $comando
fi
echo "+++ Status: OK"

# Compiling and link
echo "Compiling modules:"
echo "+++ Initializing..."
comando="gfortran -c ../src/pside.f"; eval $comando; echo -e "\t$comando"
comando="gfortran -c ../src/psidea.f"; eval $comando; echo -e "\t$comando"
comando="gfortran -c ../src/report.f"; eval $comando; echo -e "\t$comando"
comando="gfortran -c ../src/psided.f"; eval $comando; echo -e "\t$comando"
# Building executable
echo "+++ Linking..."
comando="gfortran -o ../dotest ../$1 psided.o pside.o psidea.o report.o"
eval $comando; echo -e "\t$comando"
echo "+++ Status: OK"

