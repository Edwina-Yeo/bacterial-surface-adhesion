

Vs=100
Dr=1
kappa=1
rep=1
beta=0

for ((rep=0; rep<4;rep++))
do

julia main-rate.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &

done

Vs=22
Dr=0.1
kappa=1
beta=0

for ((rep=0; rep<4;rep++))
do


julia main-rate.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &

done


Vs=70
Dr=2
kappa=1
beta=0

for ((rep=0; rep<4;rep++))
do


julia main-rate.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &

done


Vs=40
Dr=2
kappa=1
beta=0.25

for ((rep=0; rep<4;rep++))
do


julia main-rate.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &

done

Vs=40
Dr=2
kappa=1
beta=0.88

for ((rep=0; rep<4;rep++))
do


julia main-rate.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &

done

Vs=22
Dr=1
kappa=1
beta=0.1

for ((rep=0; rep<4;rep++))
do


julia main-rate.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &

done



Vs=40
Dr=4
kappa=1
beta=0.4

for ((rep=0; rep<4;rep++))
do


julia main-rate.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &

done

Vs=10
Dr=0.5
kappa=1
beta=0

for ((rep=0; rep<4;rep++))
do


julia main-rate.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &

done

