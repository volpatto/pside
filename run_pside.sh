if [ -z "$1" ]; then
    comando="time ./dotest"; eval $comando
else
    export OMP_NUM_THREADS=$1
    echo Running with number of threads: $OMP_NUM_THREADS
    comando="time ./dotest"; eval $comando
fi
