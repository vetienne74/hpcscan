
# 3D acoustic modeling in the SEG/EAGE Overthrust model

# if not present, download P-wave velocity file
VPFILE=fvp     
if [ -f $VPFILE ]; then
    echo "File $VPFILE exists."
else
    echo "File $VPFILE does not exist. Downloading it..."
    wget https://www.geoazur.fr/WIND/pub/nfs/FWI-DATA/GEOMODELS/Overthrust/fvp
    echo "Download completed."
fi

# velocity model grid size and sampling
dataPath=.
n1=187
n2=801
n3=101
d1=25
d2=25
d3=25

# clean dir
sh clean_dir.sh

# run hpcscan
mpirun -n 1 ../..//bin/hpcscan -testCase Modeling -n1 $n1 -n2 $n2 -n3 $n3 -d1 $d1 -d2 $d2 -d3 $d3 \
		     -dim 3 -ntry 1 -tmax 2.0 -snapDt 0.25 -freqMax 20.0 \
		     -modelVpFile ${dataPath}/fvp \
		     -writeGrid

# display results
sh display3DResults.sh

# end
