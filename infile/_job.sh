#! /bin/bash
# Program
i=1
while [ "$i" != "100" ]
do 
    cp in in_$i
    sed 3s/1/$i/g -i in_$i
    sed 4s/1/$i/g -i in_$i
    sed 5s/0/$i/g -i in_$i
    sed 6s/1/$(($i+1))/g -i in_$i
    i=$(($i+1))
done
