figure
clear all


set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultTextFontname', 'CMU Serif')
set(groot,'DefaultAxesFontName', 'CMU Serif')
set(groot,'DefaultLineLineWidth',2)

% set(0, 'DefaultFigureRenderer', 'painters');
set(groot,'DefaultAxesFontSize',22 ...
    );
set(0, 'DefaultFigureRenderer', 'opengl');
close all
gam=logspace(-1.5,3,20)

x1=0.1
x2=3
L=750
beta="0.88"
beta=0.88
U=gam*L
Vs=40.00;% micron/s
Vss=Vs./U
D_r=2
Pers=U./(D_r*L)%=gam/dr
D_r="2"
S=load("data/sd-Vs"+string(Vs)+"-dr-"+string(D_r)+"-beta-"+beta ...
    +".txt")'

J=load("data/J-Vs"+string(Vs)+"-dr-"+string(D_r)+"-beta-"+beta ...
    +".txt")'
ii=17
% 


figure1=figure('units','inch','position',[0,0,7.5,6 ...
    ]);
factor=@(b,per) -(128 + (-4 + b)*(-2 + b)*per.^2)./(-16 + (-4 + b^2)*per.^2) 


errorbar(Pers(1:ii),J./(Vss(1:ii).^(4/3).*factor(beta,Pers(1:ii)).^(2/3)),S./(Vss(1:ii).^(4/3).*factor(beta,Pers(1:ii)).^(2/3)),'.','Color',[0.929411764705882 0.694117647058824 0.125490196078431],HandleVisibility='off',LineWidth=2,MarkerSize=20,DisplayName='$\mathcal{V}_s=22\mu$m/s, $D_r$=1s$^{-1}$, $\beta$=0.88'); hold on

Per_vec=logspace(-3,3,100);
Pe_eff_plot=@(per) 4*(per)./((16 + per.^2))
loglog(Per_vec,1.05e7*Pe_eff_plot(Per_vec).^(2/3),'k',DisplayName='Active L\''ev$\hat{\mathrm{e}}$que theory')


set(gca,'Xscale','log') % yeah, that is indeed easy!

% set(gca,'Yscale','log') % yeah, that is indeed easy!
xlabel("Rotational P\'eclet number $Pe_r$")
ylabel('Scaled bacterial adhesion rate')%, $\int\,{J}\,\mathrm{d}x$')
title('(e)')

exportgraphics(figure1,'flux0.88.pdf')

