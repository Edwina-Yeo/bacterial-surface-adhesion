import numpy as np
import matplotlib.pyplot as plt
# from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata
from scipy.special import gamma, gammainc
import math


import numpy as np
import matplotlib.pyplot as plt

# # Creates Fig 2b.

# update pyplot settings for all plots.
plt.rcParams.update({
    "text.usetex": True,
    "font.size": 12, "font.family": "serif",
    "font.serif": ["Palatino"]})



jj=0

fig, ax = plt.subplots(1,1,figsize=(6, 4))
Perstring=['0.1','2.5','10','100']
Per=1.0
Pe_r=Per

beta=0.0
Vs = 0.01
dip = 0.0
Per=1.0
Pe_eff = 1/( Vs ** 2 * Per / (2 * (1 + (1 / 4) * Per ** 2)))

eps=Pe_eff**(-1/3)

def func(x, y,eps):
    return gammainc(1/3, (y/eps)**3 /(9* x))

# Generate values for x and y
x = np.linspace(0.01, 3, 300)
y = np.linspace(0.0, 0.5, 300)

X, Y = np.meshgrid(x, y)
Z = func(X, Y,eps)
plt.figure()

cs = plt.contour(X,Y,Z ,[0.5])
p = cs.collections[0].get_paths()[0]
v = p.vertices
x_50 = v[:,0]
y_50 = v[:,1]


# enerate values for x and yfrom mpl_toolkits.axes_grid1.inset_locator import inset_axes
x = np.linspace(0.01, 3, 300)
y = np.linspace(0.01, 0.4, 300)

X, Y = np.meshgrid(x, y)

make_plots = True
Np = 2
dt = 0.01
Np = 2
dt = 0.01
T = 300
t = 0
Ly = 0.8
Lx = 3

Nbx = 150
Nby = 200
Nt = 50
L_tol = 0.1
dt_save=200
N_arrive=100
num_in_box=2*N_arrive*Lx/(Ly*dt*Nbx*Nby) # working for Ly and Lx varying

x_edges = np.linspace(0, Lx, Nbx+1)  # Assuming some values for x_edges
y_edges = np.linspace(0,Ly, Nby+1)  # Assuming some values for y_edges
X,Y=np.meshgrid(x_edges,y_edges)

k_plot = int(T / (dt * 10))

print(k_plot)



# k = k_plot - 1
# rho = np.loadtxt(f"data/rho{k}-Vs{Vs}-beta-{beta}-Per-{Pe_r}.txt")
# total = np.ones(
#     (len(x_edges) - 1, len(y_edges) - 1))  # average counter: not all datasets will feature data in each point.

# k_av = 1
# th = np.pi * np.ones((len(x_edges) - 1, len(y_edges) - 1))

# # print(T / (dt * 10), k_plot)
# for k in range(1000, k_plot):
#     # print(k)

#     rho += np.loadtxt(f"data/rho{k}-Vs{Vs}-beta-{beta}-Per-{Pe_r}.txt" )

#     total_new = np.loadtxt(f"data/total{k}-Vs{Vs}-beta-{beta}-Per-{Pe_r}.txt")

#     for j in range(len(y_edges) - 1):
#         for i in range(len(x_edges) - 1):
#             if total_new[i, j] > 0.0:
#                 total[i, j] += total_new[i, j]

#     k_av += 1

# rho = rho / (k_av * num_in_box)

# plt.figure()
# cs = plt.contour(X[1:, 1:], Y[1:, 1:], gaussian_filter(rho, sigma=3.5), [0.5])
# #
# p = cs.collections[0].get_paths()[0]
# v = p.vertices
# x_a = v[:, 0]
# y_a = v[:, 1]



# fig2=plt.figure()
# cs = plt.contour(X, Y, func(X,Y,eps), [0.5],color='none')

# p = cs.collections[0].get_paths()[0]
# v = p.vertices
# x_f = v[:, 0]
# y_f = v[:, 1]
# plt.close(fig2)

# im0=ax.pcolormesh(X, Y, func(X,Y,eps), cmap='coolwarm',vmin=0, vmax=1.1)

# # # plt.title('Agent based sol.,  beta='+str(beta)+'a='+dip_str)
# ax.plot(x_a, y_a, 'r-',markersize=1,label='50\% agent-based')
# ax.plot(x_f,y_f,'k-',label='50\% continuum')

# #inset_ax = inset_axes(ax[0], width="40%", height=1.3, loc="upper left")
# inset_ax= ax.inset_axes([.55, .45, .4, .2])  # [x, y, width, height] w.r.t. ax
# inset_ax.pcolormesh(X[:,:], Y[:,:], rho[:,:], cmap='coolwarm',vmin=0, vmax=1.1)
# inset_ax.set_ylim([0, 0.15])

# inset_ax.plot(x_a, y_a, 'r-',markersize=1,label='50\% agent-based')
# inset_ax.plot(x_f,y_f,'k-',label='50\% continuum')
# inset_ax.set_title('Agent-based density')

# plt.xticks([]); plt.yticks([])  # strip ticks, which collide w/ main ax
# ax.set_ylim([0,0.15])
# ax.set_title('Continuum bacteria density')

# ax.set_xlabel('$x$')
# plt.colorbar(im0, ax=ax)

# if jj == 0:
#     ax.legend()
#     ax.set_ylabel('$y$')


# jj+=1



# plt.show()
# fig.savefig('flux-cont.jpg',bbox_inches='tight',dpi=500)
