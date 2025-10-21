
Vs=22
Dr=1
kappa=0.8
beta=0.1

for ((rep=0; rep<4;rep++))
do


julia main-rate.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &

done



Vs=22
Dr=1
kappa=0.6
beta=0.1

for ((rep=0; rep<4;rep++))
do


julia main-rate.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &

done



Vs=22
Dr=1
kappa=0.4
beta=0.1

for ((rep=0; rep<4;rep++))
do


julia main-rate.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &

done


