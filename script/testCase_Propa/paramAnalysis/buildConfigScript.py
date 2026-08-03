
script_template = 'templateScript.sh'

n1Range=[200,250,300,350,400,500,600,700,800,900,1000]
fdOrderRange=['2','4','6','8','10','12','14','16']
tmax='10.0'
snapDt='1.0'
dt='0.1'
nmode='150'
xmax=1000
#dim='3'
dim='2'
ntry='1'
#ratioCFLRange=['1.0','0.5','0.2','0.1']
ratioCFLRange=['1.0']
#testModeRange=['NEC_SCA']
testModeRange=['Baseline']
#propagatorRange=['Ac2Standard','Ac2SplitComp']
propagatorRange=['Ac2Standard']

#-----------------------------------------------------------------
configName='config'
ii = 1
for propagator in propagatorRange:
    for testMode in testModeRange:
        for n1 in n1Range:
            for fdOrder in fdOrderRange:
                for ratioCFL in ratioCFLRange:

                    # create new script
                    new_name = configName+str(ii)
                    script_new = script_template.replace('template', new_name)                

                    # skip some unwanted config
                    #if propagator == 'Ac2Standard' and testMode == 'NEC_SCA':
                    #    continue

                    # open files
                    f2 = open(script_new, 'w')
                    f1 = open(script_template, 'r')

                    # compute spatial sampling
                    d1 = xmax / (n1-1)

                    # create new file
                    for line in f1:
                        new_line = line.replace('_n1_', str(n1))
                        new_line = new_line.replace('_n2_', str(n1))
                        new_line = new_line.replace('_n3_', str(n1))
                        new_line = new_line.replace('_d1_', str(d1))
                        new_line = new_line.replace('_d2_', str(d1))
                        new_line = new_line.replace('_d3_', str(d1))
                        new_line = new_line.replace('_tmax_', tmax)
                        new_line = new_line.replace('_ntry_', ntry)
                        new_line = new_line.replace('_fdOrder_', fdOrder)
                        new_line = new_line.replace('_snapDt_', snapDt)
                        new_line = new_line.replace('_testMode_', testMode)
                        new_line = new_line.replace('_propagator_', propagator)
                        new_line = new_line.replace('_nmode_', nmode)
                        new_line = new_line.replace('_dt_', dt)
                        new_line = new_line.replace('_ratioCFL_', ratioCFL)
                        new_line = new_line.replace('_dim_', dim)

                        f2.write(new_line)

                    # close files
                    f1.close()
                    f2.close()

                    ii = ii + 1

# END
