
N_arrive=100
Ly=1.5
dt=0.01
Lx=3
cin=N_arrive*2/(dt*(Ly)^2) % converts dimensional density to number 
% cin=2*N_arrive*Lx/(dt*Ly) % using the total number in the box 
T_window=2500 % size of time window.
x_int=2.17518 % integral of x^(-1/3) between 0.5 and 3

coef=x_int*3^(1/3)*cin*T_window/gamma(1/3)





