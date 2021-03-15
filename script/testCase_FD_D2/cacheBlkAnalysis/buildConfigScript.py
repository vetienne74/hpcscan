
script_template = 'templateScript.sh'

n1Range=['1000']
fdOrderRange=['4','8']
dim='1'
ntry='5'
cb1Range=[9999]
# from 1 to 10
cb2Range=range(1,11,1)
cb3Range=range(1,11,1)

#-----------------------------------------------------------------
configName='config'
ii = 1
for n1 in n1Range:
    for fdOrder in fdOrderRange:
        for cb1 in cb1Range:
            for cb2 in cb2Range:
                for cb3 in cb3Range:

                    # create new script
                    new_name = configName+str(ii)
                    script_new = script_template.replace('template', new_name)                

                    # skip some unwanted config
                    # none at present
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
                        new_line = new_line.replace('_ntry_', ntry)
                        new_line = new_line.replace('_fdOrder_', fdOrder)
                        new_line = new_line.replace('_dim_', dim)
                        new_line = new_line.replace('_cb1_', str(cb1))
                        new_line = new_line.replace('_cb2_', str(cb2))
                        new_line = new_line.replace('_cb3_', str(cb3))

                        f2.write(new_line)

                    # close files
                    f1.close()
                    f2.close()

                    ii = ii + 1

# END
