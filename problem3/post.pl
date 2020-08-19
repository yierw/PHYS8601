#!/usr/bin/perl

#----------------------------
# It is a perl script to calculate the mean 
# and error bar of the simulation results
#----------------------------
$file="runs.dat";
$col=4;
$row=10;

for($i=0;$i<$col;$i++)
{
    $sum[$i]=0;
    $sum2[$i]=0;
}

open(FD,$file);
@sim=<FD>;
close(FD);

foreach(@sim)
{
    @sep=split(' ',$_);
    $i=0;
    foreach(@sep)
    {
        $sum[$i]+=$_;
        $sum2[$i]+=$_*$_;
        $i++;
   }
}


foreach(@sum){push(@mean,$_/$row);}

for ($i=0;$i<$col;$i++)
{
    $s=sqrt(($sum2[$i]-$row*@mean[$i]**2)/($row-1.0));
    push(@sig,$s);
}

for ($i=0;$i<$col;$i++)
{
   printf "%10.6f %10.6f  ",@mean[$i],@sig[$i];
}
printf "\n "
