
cp ../../src/grid_Cuda.cu .

hipconvertinplace-perl.sh grid_Cuda.cu

sed -i -e 's/Cuda/Hip/g' grid_Cuda.cu

sed -i -e 's/CUDA/HIP/g' grid_Cuda.cu

mv grid_Cuda.cu ../../src/grid_Hip.hip

rm grid_Cuda.cu.prehip

