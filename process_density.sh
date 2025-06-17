#!/bin/bash



#Pe=24000
#beta=0.
#Vs=0.0
#dip=0.0
#for Per in 1; do
#  julia main_passive.jl  "${Pe}" "${beta}" "${Vs}" "${dip}" "${Per}"&
#done

#Pe=2000
#beta=0.0
#dip=0.0
#Per=.1
#for Vs in 0.05; do
# # julia main.jl  "${Pe}" "${beta}" "${Vs}" "${dip}" "${Per}"&
#done
#Pe=2000
#beta=0.5
#dip=0.0
#Per=1.0
#for Vs in 0.05; do
#  julia main.jl  "${Pe}" "${beta}" "${Vs}" "${dip}" "${Per}"&
#done



#Pe=2000
#beta=0.0

#done
Pe=2000
beta=0
dip=0.0
#!/bin/bash

# Arrays with normal spacing between elements
Pers=(0.0316227766016838 0.0545559478116852 0.0941204967268067 0.162377673918872 0.280135676119887 0.483293023857175 0.833782223471789 1.43844988828766 2.48162892283683 4.28133239871940 7.38619982207936 12.7427498570313 21.9839264886229 37.9269019073225 65.4318912971297 112.883789168469 194.748303990876 335.981828628378 579.639395338497 1000)

Vss=(0.927601446982725 0.537674341844181 0.311657230395585 0.180648808579368 0.104711166173565 0.0606947170460338 0.0351810490888042 0.0203923220212084 0.0118201932059212 0.00685144964266436 0.00397137012806613 0.00230196257969762 0.00133430819778777 0.000773417596960905 0.000448303308246572 0.000259854258520291 0.000150621765284835 8.73063089545187e-05 5.06061761316331e-05 2.93333333333333e-05)


length=20

for ((j=0; j<20;j++))
do 
	Vs=${Vss[$j]}
	Per=${Pers[$j]}

	julia main.jl  "${Pe}" "${beta}" "${Vs}" "${dip}" "${Per}" &
	if ((j==9)); then 
		echo "waiting to finish $beta $j+1 jobs"
		wait
		fi 	

done
wait 	

echo "now solving for ecoli"
beta=0.88
for ((j=0; j<20;j++))
do 
	Vs=${Vss[$j]}
	Per=${Pers[$j]}

	julia main.jl  "${Pe}" "${beta}" "${Vs}" "${dip}" "${Per}" &
	if ((j==9)); then 
		echo "waiting to finish $beta $j+1 jobs"
		wait
		fi 	

done
wait 	

		




