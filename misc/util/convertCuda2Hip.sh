
# convert grid_Cuda.cu to grid_Hip.hip
cp ../../src/grid_Cuda.cu .
hipconvertinplace-perl.sh grid_Cuda.cu
sed -i -e 's/Cuda/Hip/g' grid_Cuda.cu
sed -i -e 's/CUDA/HIP/g' grid_Cuda.cu
mv grid_Cuda.cu ../../src/grid_Hip.hip
rm grid_Cuda.cu.prehip

# convert grid_Cuda.h to grid_Hip.h
cp ../../src/grid_Cuda.h .
sed -i -e 's/Cuda/Hip/g' grid_Cuda.h
sed -i -e 's/CUDA/HIP/g' grid_Cuda.h
mv grid_Cuda.h ../../src/grid_Hip.h


