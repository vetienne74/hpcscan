
n1=187
n2=801
d1=25
d2=25

# display input velocity model
ximage < ./fvp n1=$n1 n2=$n2 d1=$d1 d2=$d2 \
       title='Overthrust Vp model' \
       legend=1 hbox=$n1 wbox=$n2 &

# display computed snapshots
#ximage < ModelingForwardPrn.global.grid.bin n1=187 \
#       title='Computed snapshots' \
#       perc=99 &

# display computed traces
ximage < ModelingForwardPrn.trace.bin n1=873 \
       title='Computed traces' \
       perc=99 &
