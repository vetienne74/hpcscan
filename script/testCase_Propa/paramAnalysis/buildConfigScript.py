
script_template = 'templateScript.sh'

n1Range=['500','600','700','800','900','1000']
fdOrderRange=['4','6','8','10','12']
tmax='1.0'
snapDt='0.1'
dt='0.0'
nmode='130'
#dim='3'
dim='1'
ntry='1'
ratioCFLRange=['1.0','0.5','0.2','0.1']
testModeRange=['NEC_SCA']
propagatorRange=['Ac2Standard','Ac2SplitComp']

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

                    # create new file
                    for line in f1:
                        new_line = line.replace('_n1_', n1)
                        new_line = new_line.replace('_n2_', n1)
                        new_line = new_line.replace('_n3_', n1)
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
