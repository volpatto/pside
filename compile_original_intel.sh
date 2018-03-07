###############################################################
#                                                             #
# Brief description:                                          #
# This script compile and link PSIDE source codes with ifort  #
# compiler. User can employ 2 using modes:                    # 
#   (i) Auto-parallel region: In this case, user provides as  #
#       second arg the numbers of OMP threads;                #
#   (ii) Regular: ifort standard compliling and linking. No   #
#        additional args                                      #
# Summarizing, the use cases is:                              #
#   ./compile_intel.sh test_problem Nt                        #
# where test_problem is the fortran source file which descri- #
# bes the problem to be solved by PSIDE and Nt means the num- #
# of OMP threads, if user wants to employ auto-parallel mode. #
# This last arg can be ommited.                               #
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

# Verifing ifort version
echo "Compiler version check:"
echo -e "+++ Version:\n"
comando="ifort -V"; eval $comando
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
if [ -z "$2" ]; then
    # Compiling modules
    echo "+++ Initializing..."
    comando="ifort -c ../original/pside.f"; eval $comando; echo -e "\t$comando"
    comando="ifort -c ../original/psidea.f"; eval $comando; echo -e "\t$comando"
    comando="ifort -c ../original/report.f"; eval $comando; echo -e "\t$comando"
    comando="ifort -c ../original/psided.f"; eval $comando; echo -e "\t$comando"
    # Building executable
    echo "+++ Linking..."
    comando="ifort -o ../dotest ../$1 psided.o pside.o psidea.o report.o"
    eval $comando; echo -e "\t$comando"
    echo "+++ Status: OK"
else
    export OMP_NUM_THREADS=$2
    echo +++ Number of threads: $OMP_NUM_THREADS
    # Compiling modules
    echo "+++ Initializing..."
    comando="ifort -c -parallel ../original/pside.f"; eval $comando; echo -e "\t$comando"
    comando="ifort -c -parallel ../original/psidea.f"; eval $comando; echo -e "\t$comando"
    comando="ifort -c -parallel ../original/report.f"; eval $comando; echo -e "\t$comando"
    comando="ifort -c -parallel ../original/psided.f"; eval $comando; echo -e "\t$comando"
    # Building executable
    echo "+++ Linking..."
    comando="ifort -parallel -o ../dotest ../$1 psided.o pside.o psidea.o report.o"
    eval $comando; echo -e "\t$comando"
    echo "+++ Status: OK"
fi

