

set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultTextFontname', 'CMU Serif')
set(groot,'DefaultAxesFontName', 'CMU Serif')
set(groot,'DefaultLineLineWidth',2)

set(groot,'DefaultAxesFontSize',22);
set(0, 'DefaultFigureRenderer', 'opengl');
close all
figure1=figure('units','inch','position',[0,0,6,6 ...
    ]);

gam=logspace(-1.5,3,20)
ii=16



% Calculate the scaling from agent based simulations 
N_arrive=100;
Ly=1.5;
dt=0.01;
Lx=3;
cin=N_arrive*2/(dt*(Ly)^2); % converts dimensional density to number 
% cin=2*N_arrive*Lx/(dt*Ly) % using the total number in the box 
T_window=2500; % size of time window.
x_int=2.17518; % integral of x^(-1/3) between 0.5 and 3
coeff_line=x_int*3^(1/3)/gamma(1/3);
coeff_agent=cin.*T_window;


factor=@(b,per) (32 + (2- b)*(1- b)*per.^2)./(16 + (4 - b^2)*per.^2) 
Pe_eff_plot=@(per) 4*(per)./((16 + per.^2))
Pe_isotropic_plot=@(per) (per)/4

L=750
U=gam*L

gam=logspace(-1.5,3,20)


beta_str="0.88"
Vs=40.00;% micron/s
D_r=2
D_r_str="2"
beta=0.88
colour='r'
add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff_agent)

Per_vec=logspace(-3,3,100);

loglog(Per_vec,       coeff_line*Pe_eff_plot(Per_vec).^(2/3),'k',DisplayName='Active L\''ev$\hat{\mathrm{e}}$que theory')

set(gca,'Xscale','log') % yeah, that is indeed easy!
xlabel("Rotational P\'eclet number $Pe_r$")
ylabel('Scaled total bacterial adhesion rate, $\bar{J}$')%, $\int\,\hat{J}/V_s^{4/3}\,\mathrm{d}x$')
legend(location='northwest')

xlim([1e-3,1e3])
ylim([5e-3,2.3])




exportgraphics(figure1,'flux'+beta_str+'.pdf')


function add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff)
factor=@(b,per) (32 + (2- b)*(1- b)*per.^2)./(16 + (4 - b^2)*per.^2) ;
ii=16;

Vss=Vs./U(1:ii);
Pers=U(1:ii)./(D_r*L);%=gam/dr
S=load("data/sd-Vs"+string(Vs)+"-dr-"+string(D_r)+"-beta-"+beta_str ...
    +".txt")';

J=load("data/J-Vs"+string(Vs)+"-dr-"+string(D_r)+"-beta-"+beta_str ...
    +".txt")';

JJ=J(1:ii).*gam(1:ii)/2500;
Err=S(1:ii).*gam(1:ii)/2500; % scaled by the time window we are av

errorbar(Pers(1:ii),J(1:ii)./(coeff*Vss(1:ii).^(4/3).*factor(beta,Pers(1:ii)).^(2/3)), ...
    S(1:ii)./(coeff*Vss(1:ii).^(4/3).*factor(beta,Pers(1:ii)).^(2/3)),'.', ...
    'Color',colour,HandleVisibility='on', ...
    LineWidth=2,MarkerSize=20,DisplayName='Agent-based data, $\beta$='+string(beta)); hold on

end


