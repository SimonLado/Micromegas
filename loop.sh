var2_values=(0.05 0.10 0.12 0.14 0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30)

for var2 in ${var2_values[@]}; do
    a=1
    while [ $a -le 500 ]; do
        let c=`qstat -u ladogana | grep ladogana | wc -l`
        echo "jobs= $c"
        if (( c < 200 )); then
            qsub -v "var1=480,var2=${var2},var3=${a}" mmjob.sh
            let a=a+1
            echo $a
            sleep .5
        else
            echo "Waiting........"
            sleep 10
        fi
    done
done