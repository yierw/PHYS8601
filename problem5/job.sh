name="reweightingcv.dat"
a=0
while [ "$a" -lt 10 ]    # 10 independent runs
do
   echo $a
   ./run.x
   sed '2d' $name> $a.dat
   a=`expr $a + 1`
   sleep 1
done





