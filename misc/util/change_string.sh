
dir='../../script'

echo 'change_string start!'

for file in $(find ${dir} -name '*')
do
    echo 'process ' ${file}
    sed -i -e 's/KMP_AFFINITY=scatter,1,0,granularity=fine/KMP_AFFINITY=granularity=fine,compact/g' $file
done

echo 'change_string done!'

exit

#-----------------------------------------------------------------------------------

# change string in sources files
#dir='../src'
#dir='../script'
dir='../hackathonTestCases/testCase_Propa/'
pattern='*'

echo 'change_string start!'

for file in $(ls ${dir}/${pattern}):
do

    echo 'process ' ${file}   
    sed -i -e 's/propaName/propagator/g' $file
done

echo 'change_string done!'

exit

#-----------------------------------------------------------------------------------

# change string in sources files
dir='../src'
#dir='../include'
pattern='*'

echo 'change_string start!'

for file in $(ls ${dir}/${pattern}):
do

    echo 'process ' ${file}   
    sed -i -e 's/printDebug(Display_type, /printDebug(/g' $file
done

echo 'change_string done!'

exit
