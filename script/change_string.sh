
#-----------------------------------------------------------------------------------

# change string in sources files
dir='../src'
#dir='../script'
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
