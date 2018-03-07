###############################################################
#                                                             #
# Brief description:                                          #
# This script compile and link PSIDE source codes with PGI    #
# compiler. User can employ 2 using modes:                    # 
#   (i) Code optimization: In this case, user sign with an    #
#       arg and compiler optimizes the code;                  #
#   (ii) Regular: pgfortran standard compliling and linking.  #
#        No additional args                                   #
# Summarizing, the use cases is:                              #
#   ./compile_pgi.sh test_problem flag                        #
# where test_problem is the fortran source file which descri- #
# bes the problem to be solved by PSIDE and flag means the    #
# sign of optimization.                                       #
# This last arg can be ommited. Sign should be number 1.      #
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

# Verifing pgfortran version
echo "Compiler version check:"
echo -e "+++ Version:"
comando="pgfortran -V"; eval $comando
if [ $? -ne 0 ]; then
    echo "+++ Compiler can not be checked"
    echo "+++ Status: Failed"
    exit 1
fi
echo -e "\n+++ Status: OK\n"

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
if [ -z "$2" ]; then
    # Compiling modules
    echo "+++ Initializing..."
    comando="pgfortran -c ../original/pside.f"; eval $comando; echo -e "\t$comando"
    comando="pgfortran -c ../original/psidea.f"; eval $comando; echo -e "\t$comando"
    comando="pgfortran -c ../original/report.f"; eval $comando; echo -e "\t$comando"
    comando="pgfortran -c ../original/psided.f"; eval $comando; echo -e "\t$comando"
    # Building executable
    echo "+++ Linking..."
    comando="pgfortran -o ../dotest ../$1 psided.o pside.o psidea.o report.o"
    eval $comando; echo -e "\t$comando"
    echo "+++ Status: OK"
elif [ $2 -eq 1 ]; then
    # Compiling modules
    echo "+++ Initializing..."
    comando="pgfortran -c -fast ../original/pside.f"; eval $comando; echo -e "\t$comando"
    comando="pgfortran -c -fast ../original/psidea.f"; eval $comando; echo -e "\t$comando"
    comando="pgfortran -c -fast ../original/report.f"; eval $comando; echo -e "\t$comando"
    comando="pgfortran -c -fast ../original/psided.f"; eval $comando; echo -e "\t$comando"
    # Building executable
    echo "+++ Linking..."
    comando="pgfortran -fast -o ../dotest ../$1 psided.o pside.o psidea.o report.o"
    eval $comando; echo -e "\t$comando"
    echo "+++ Status: OK"
else
    echo "Error: The compiling mode was not specified"
    exit 1
fi

