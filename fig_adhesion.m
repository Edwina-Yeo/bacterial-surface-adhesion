figure
clear all

set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultTextFontname', 'CMU Serif')
set(groot,'DefaultAxesFontName', 'CMU Serif')
set(groot,'DefaultLineLineWidth',2)

set(groot,'DefaultAxesFontSize',19);
set(0, 'DefaultFigureRenderer', 'opengl');
close all
figure1=figure('units','inch','position',[0,0,18,6 ...
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
coeff=x_int*3^(1/3)*cin*T_window/gamma(1/3);



factor=@(b,per) (32 + (2- b)*(1- b)*per.^2)./(16 + (4 - b^2)*per.^2) 
Pe_eff_plot=@(per) 4*(per)./((16 + per.^2))
Pe_isotropic_plot=@(per) (per)/4

L=750
U=gam*L

gam=logspace(-1.5,3,20)


beta_str="0.0"
Vs=100.00;% micron/s
D_r=1
D_r_str="1"
beta=0
colour=[1 0.5 0.95]
add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff)


beta_str="0.0"
Vs=22.00;% micron/s
D_r=0.1
D_r_str="0.1"
beta=0
colour=[0.0745098039215686 0.623529411764706 1]
add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff)


beta_str="0.0"
Vs=70.00;% micron/s
D_r=2
D_r_str="2"
beta=0
colour=[0.392156862745098 0.831372549019608 0.0745098039215686]
add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff)


beta_str="0.25"
Vs=40.00;% micron/s
D_r=2
D_r_str="2"
beta=0.25
colour=[0.929411764705882 0.694117647058824 0.125490196078431]
add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff)


beta_str="0.1"
Vs=22.00;% micron/s
D_r=1
D_r_str="1.0"
beta=0.1
colour='r'
add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff)

beta_str="0.5"
beta=0.5
U=gam*L
Vs=40.00;% micron/s
Vss=Vs./U
D_r=4
Pers=U./(D_r*L)%=gam/dr
D_r_str="4"
colour=[0 0.447058823529412 0.741176470588235]
add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff)


beta_str="0.0"
U=gam*L
Vs=10.00;% micron/s
D_r=0.5
D_r_str="0.5"
beta=0
colour=[1 0.411764705882353 0.16078431372549]
add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff)

Per_vec=logspace(-3,3,100);

loglog(Per_vec,       Pe_eff_plot(Per_vec).^(2/3),'k',DisplayName='Active L\''ev$\hat{\mathrm{e}}$que theory')
loglog(Per_vec, Pe_isotropic_plot(Per_vec).^(2/3),'-.','Color',[0.5,0.5,0.5],LineWidth=3,DisplayName='Isotropic theory')
loglog(1,1,'-.','Color',[1,1,1],LineWidth=3,DisplayName='Agent-based data')

set(gca,'Xscale','log') % yeah, that is indeed easy!
xlabel("Rotational P\'eclet number $Pe_r$")
ylabel('Scaled total bacterial adhesion rate')%, $\int\,\hat{J}/V_s^{4/3}\,\mathrm{d}x$')
legend(location='northwest')

xlim([5e-3,1e3])
ylim([5e-3,1])

subplot(1,2,1)
set(gca,'Xscale','log') % yeah, that is indeed easy!
set(gca,'Yscale','log') % yeah, that is indeed easy!
xlabel('Fluid shear rate $\dot{\gamma} $ (s$^{-1}$)')
ylabel('Total bacterial adhesion rate (s$^{-1}$)')%, $\int\,{J}\,\mathrm{d}x$')

legend(location='southwest')


loglog(gam(1:6),3e7*gam(1:6).^(1/3),'k-.',HandleVisibility='off')
loglog(gam(11:16),4.5e9*gam(11:16).^(-1 ...
    ),'k-.',HandleVisibility='off')
%
ylim([8e4,9e8])
 xlim([2e-2,gam(ii+1)])
exportgraphics(figure1,'flux.pdf')


function add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff)
factor=@(b,per) (32 + (2- b)*(1- b)*per.^2)./(16 + (4 - b^2)*per.^2) ;
ii=16;

Vss=Vs./U(1:ii);
Pers=U(1:ii)./(D_r*L);%=gam/dr
S=load("data/sd-Vs"+string(Vs)+"-dr-"+string(D_r)+"-beta-"+beta_str ...
    +".txt")';

J=load("data/J-Vs"+string(Vs)+"-dr-"+string(D_r)+"-beta-"+beta_str ...
    +".txt")';

subplot(1,2,1)
JJ=J(1:ii).*U(1:ii);
Err=S(1:ii).*U(1:ii);
errorbar(gam(1:ii),JJ,Err,'.','Color',colour, ...
    LineWidth=2,MarkerSize=20,DisplayName='$\mathcal{V}_s='+string(Vs)+'\mu$m/s, $D_r$='+D_r_str...
    +'s$^{-1}$, $\beta$='+string(beta)); hold on
subplot(1,2,2)
errorbar(Pers(1:ii),J(1:ii)./(coeff*Vss(1:ii).^(4/3).*factor(beta,Pers(1:ii)).^(2/3)), ...
    S(1:ii)./(coeff*Vss(1:ii).^(4/3).*factor(beta,Pers(1:ii)).^(2/3)),'.', ...
    'Color',colour,HandleVisibility='off', ...
    LineWidth=2,MarkerSize=20); hold on

end
