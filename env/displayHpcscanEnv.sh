
echo hpcscan set for `hostname` with:

echo - C++ Compiler
$HPCSCAN_CPP --version

if [ -f $HPCSCAN_CUDA ]
then
    echo - No CUDA compiler set
else
    echo - CUDA Compiler
    $HPCSCAN_CUDA --version
fi

if [ -f $HPCSCAN_HIP ]
then
    echo - No HIP compiler set
else
    echo - HIP Compiler
    $HPCSCAN_HIP --version
fi

echo - MPI library
$HPCSCAN_MPI_INVOKER --version

# custom mode
if [ -f $HPCSCAN_CUSTOM_PATH ]
then
    echo - Mode Custom disabled
else
    echo - Mode Custom enabled
fi

echo - Number of OpenMP threads $HPCSCAN_NTHREADS

echo Ready to go!

