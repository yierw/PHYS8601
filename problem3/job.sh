#!/bin/sh
echo >runs.dat
a=0
while [ "$a" -lt 10 ]    # 10 independent runs
do
   echo $a
   a=`expr $a + 1`
   ./run.x>>runs.dat
   sleep 1  
done

./post.pl


