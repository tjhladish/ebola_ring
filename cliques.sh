for i in `seq 1 100`;
do
    ./ebola_net abc_ebola_rings-hh_cliques.json --simulate -n 20 2>> ./output/part{$1}.err >> ./output/part{$1}.out
done
