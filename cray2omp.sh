# Check if input file is supplied
if [ -z "$1" ]; then
    echo "Error: Input file was not supplied"
    exit 1
fi

# Replacing DOALL to equiv PARALLEL DO of Open MP
comando="sed -ir 's/DOALL/PARALLEL DO/' $1"; eval $comando
# Replacing CMIC$+ to equiv C$OMP& of Open MP
comando="sed -ir 's/CMIC\$\+/c\$omp\&/' $1"; eval $comando
# Replacing CMIC$ to equiv C$OMP of Open MP
comando="sed -ir 's/CMIC\$/c\$omp/' $1"; eval $comando
