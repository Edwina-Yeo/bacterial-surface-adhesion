


 # Processing script for SI plot Fig. S4b
Vs=22
Dr=1
kappa=1
beta=0.1

for ((rep=0; rep<4;rep++))
do
julia main-rate-dip.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &
done

Vs=70
Dr=2
kappa=1
beta=0

for ((rep=0; rep<4;rep++))
do
julia main-rate-dip.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &
done



Vs=10
Dr=0.5
kappa=1
beta=0.1

for ((rep=0; rep<4;rep++))
do
julia main-rate-dip.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &
done


# Processing for SI plot Fig. S4a 

Vs=22
Dr=1
beta=0.1


kappa=0.6
for ((rep=0; rep<4;rep++))
do
julia main-rate-kappa.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &
done

kappa=0.4
for ((rep=0; rep<4;rep++))
do
julia main-rate-kappa.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &
done

kappa=1
for ((rep=0; rep<4;rep++))
do
julia main-rate-kappa.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &
done

kappa=0.8
for ((rep=0; rep<4;rep++))
do
julia main-rate-kappa.jl  "${beta}" "${Vs}"  "${Dr}" "${kappa}" "${rep}" &
done






