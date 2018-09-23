regex=".*/([0-9]+_gauss).csv"

for f in `ls ../posterior_networks/*gauss.csv`;
do
    if [[ $f =~ $regex ]]
    then
        ./layout_network $f > ${BASH_REMATCH[1]}.layout
        #echo ${BASH_REMATCH[1]}
    fi
    #echo $i
done
