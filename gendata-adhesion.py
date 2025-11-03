import numpy as np


def calc_flux(Vss,Pers,beta,Vs_val,Dr,gam):

    J=np.zeros(17,)
    sd=np.zeros(17,)
    kappa=0.1
    rep_num=0.0
    loc='data/'
    for i in range(0,17):

        t=np.linspace( 0, 2800,2)
        x=np.linspace(0.0,3,10)
        Per=Pers[i]
 
        Vs=Vss[i]
        dip=dips[i]
        xb = np.loadtxt(loc+'Vs' + str(Vs) + '-beta-' + str(beta) + '-Per-' + str(Per) +"-kappa-"+str(kappa)+"-rep-"+str(rep_num)+"-xbs.txt")
        tb = np.loadtxt(loc+'Vs' + str(Vs) + '-beta-' + str(beta)+ '-Per-' + str(Per)  +"-kappa-"+str(kappa)+"-rep-"+str(rep_num)+'-tbs.txt')

        print(np.shape(tb),np.shape(x))
        H, xedges, yedges = np.histogram2d(xb, tb, bins=(x, t))
        H = H.T
        print(H)
        ad_t=np.sum(H,axis=1)
    
        J[i]=np.mean(ad_t)

#    np.savetxt('data/J-Vs'+str(Vs_val)+'-dr-'+str(Dr)+'-beta-'+str(beta)++"-kappa-"+str(kappa)+'.txt',J)

    print('the dimensionless flux is',J)
    print('the dimensional flux is',J*gam[:17])




L=750
gam=np.logspace(-1.5,2,20)

Vs=22
Dr=1
beta=0.1

Vss = [0.927601446982725, 0.537674341844181, 0.311657230395585, 0.180648808579368, 0.104711166173565, 0.0606947170460338, 0.0351810490888042, 0.0203923220212084, 0.0118201932059212, 0.00685144964266435, 0.00397137012806613, 0.00230196257969762, 0.00133430819778777, 0.000773417596960905, 0.000448303308246572, 0.000259854258520291, 0.000150621765284835, 8.73063089545187e-05, 5.06061761316331e-05, 2.93333333333333e-05]
Pers = [0.0316227766016838, 0.0545559478116852, 0.0941204967268067, 0.162377673918872, 0.280135676119887, 0.483293023857175, 0.833782223471789, 1.43844988828766, 2.48162892283682, 4.2813323987194, 7.38619982207936, 12.7427498570313, 21.9839264886229, 37.9269019073225, 65.4318912971297, 112.883789168469, 194.748303990876, 335.981828628378, 579.639395338497, 1000]


calc_flux(Vss,Pers,beta,Vs,Dr,gam)

