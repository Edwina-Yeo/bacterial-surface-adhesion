clear all]
% Generates Fig1b and Fig2c. Figures are combined using SVG files for paper
% versions.

% Set plotting macros
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

gam=logspace(-1.5,3,20); % Shear rates 

% Calculate the scaling from agent based simulations to convert number to
% denisty
N_arrive=100;
Ly=1.5;
dt=0.01;
Lx=3;
cin=N_arrive*2/(dt*(Ly)^2); % converts dimensional density to number 
T_window=2500; % size of time window.
x_int=2.17518; % integral of x^(-1/3) between 0.5 and 3
coeff_line=x_int*3^(1/3)/gamma(1/3);
coeff_agent=cin.*T_window;

% Leveque theory for plotting. 
Pe_eff_plot=@(per) 4*(per)./((16 + per.^2))
Pe_isotropic_plot=@(per) (per)/4

L=750 % lengthscale for re-dimensionalising 
U=gam*L

gam=logspace(-1.5,3,20)


% Load agent based datasets--------------

beta_str="0.0";
Vs=100.00;
D_r=1;
D_r_str="1";
beta=0;
colour=[1 0.5 0.95];
add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff_agent)


beta_str="0.0";
Vs=22.00;
D_r=0.1;
D_r_str="0.1";
beta=0;
colour=[0.0745098039215686 0.623529411764706 1];
add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff_agent)


beta_str="0.0";
Vs=70.00;
D_r=2;
D_r_str="2";
beta=0;
colour=[0.392156862745098 0.831372549019608 0.0745098039215686];
add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff_agent)

beta_str="0.25";
Vs=40.00;
D_r=2;
D_r_str="2";
beta=0.25;
colour=[0.929411764705882 0.694117647058824 0.125490196078431];
add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff_agent)

beta_str="0.1";
Vs=22.00;
D_r=1;
D_r_str="1.0";
beta=0.1;
colour='r';
add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff_agent)


beta_str="0.4";
beta=0.4;
U=gam*L;
Vs=40.00;;
Vss=Vs./U;
D_r=4;
Pers=U./(D_r*L);
D_r_str="4";
colour=[0 0.447058823529412 0.741176470588235]
add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff_agent)


beta_str="0.0";
U=gam*L;
Vs=10.00;
D_r=0.5;
D_r_str="0.5";
beta=0;
colour=[1 0.411764705882353 0.16078431372549];
add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff_agent)


% Add Leveque theory lines over-top
Per_vec=logspace(-3,3,100);
loglog(Per_vec,       coeff_line*Pe_eff_plot(Per_vec).^(2/3),'k',DisplayName='Active L\''ev$\hat{\mathrm{e}}$que theory')
loglog(Per_vec, coeff_line*Pe_isotropic_plot(Per_vec).^(2/3),'-.','Color',[0.5,0.5,0.5],LineWidth=3,DisplayName='Isotropic theory')
loglog(1,1,'-.','Color',[1,1,1],LineWidth=3,DisplayName='Agent-based data')

set(gca,'Xscale','log') 
xlabel("Rotational P\'eclet number $Pe_r$")
ylabel('Scaled total bacterial adhesion rate, $\bar{J}$')
legend(location='northwest')
xlim([1e-3,1e3])
ylim([5e-3,1])

subplot(1,2,1)
set(gca,'Xscale','log') 
set(gca,'Yscale','log')
xlabel('Fluid shear rate $\dot{\gamma} $ (s$^{-1}$)')
ylabel('Total bacterial adhesion rate (cells/s$^{-1}$)')
legend(location='southwest')
% Add scaling lines on top
loglog(gam(1:6),3e7*gam(1:6).^(1/3)./(2500*750),'k-.',HandleVisibility='off')
loglog(gam(11:16),4.5e9*gam(11:16).^(-1 ...
    )./(2500*750),'k-.',HandleVisibility='off')
ylim([8e4,9e8]./(2500*750))
 xlim([2e-2,gam(ii+1)])
exportgraphics(figure1,'flux.pdf')


function add_agent_data(L,gam,beta_str,beta,Vs,U,D_r,D_r_str,colour,coeff)

% Scale factor S(\beta,Per,Vss) eq 18.
factor=@(b,per,Vss) (32 + (2- b)*(1- b)*per.^2)./(16 + (4 - b^2)*per.^2).*Vss.^2 ;
ii=16; % Number of points plotted 

% Create dimensionless param vectors
Vss=Vs./U(1:ii);
Pers=U(1:ii)./(D_r*L);

% Load agent based data
S=load("data/sd-Vs"+string(Vs)+"-dr-"+string(D_r)+"-beta-"+beta_str ...
    +".txt")';

J=load("data/J-Vs"+string(Vs)+"-dr-"+string(D_r)+"-beta-"+beta_str ...
    +".txt")';

subplot(1,2,1)
% Dimensional plot
JJ=J(1:ii).*gam(1:ii)/2500; % scaled by the time window we are averaging over
Err=S(1:ii).*gam(1:ii)/2500; % scaled by the time window we are averaging over
errorbar(gam(1:ii),JJ,Err,'.','Color',colour, ...
    LineWidth=2,MarkerSize=20,DisplayName='$\mathcal{V}_s='+string(Vs)+'\mu$m/s, $D_r$='+D_r_str...
    +'s$^{-1}$, $\beta$='+string(beta)); hold on
subplot(1,2,2)
% ND plot
errorbar(Pers(1:ii),J(1:ii)./(coeff.*factor(beta,Pers(1:ii),Vss(1:ii)).^(2/3)), ...
    S(1:ii)./(coeff*factor(beta,Pers(1:ii),Vss(1:ii)).^(2/3)),'.', ...
    'Color',colour,HandleVisibility='off', ...
    LineWidth=2,MarkerSize=20); hold on

end
