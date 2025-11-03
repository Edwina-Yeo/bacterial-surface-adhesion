%% Fig S4a: imperfect adhesion results 

figure1=figure('units','inch','position',[0,0,7,6 ...
    ]);

Vs=22
D_r_str='1'
beta=0.1
gam=logspace(-1.5,2,20); % Shear rates 
colour(1,:)=[1,0,0];
colour(2,:)=[0.392156862745098 0.831372549019608 0.0745098039215686];
colour(3,:)=[0.929411764705882 0.694117647058824 0.125490196078431];
colour(4,:)=[0 0.447058823529412 0.741176470588235]

i=1
J=[];
for kappa=[1.0]
J(1,:)=load("data/Vs22.0-beta-0.1-D_r-1.0-kappa-1.0-rep-0.0-rate.txt");
J(2,:)=load("data/Vs22.0-beta-0.1-D_r-1.0-kappa-1.0-rep-1.0-rate.txt");
J(3,:)=load("data/Vs22.0-beta-0.1-D_r-1.0-kappa-1.0-rep-2.0-rate.txt");
J(4,:)=load("data/Vs22.0-beta-0.1-D_r-1.0-kappa-1.0-rep-3.0-rate.txt");

errorbar(gam,mean(J,1),std(J,1),'.','Color',colour(i,:),LineWidth=2,MarkerSize=20,DisplayName='$\mathcal{V}_s='+string(Vs)+'\mu$m/s, $D_r$='+D_r_str...
    +'s$^{-1}$, $\beta$='+string(beta)+', $\kappa=$'+string(kappa));hold on
i=i+1
end

for kappa=[0.8:-0.2:0.4]
    J=[];
J(1,:)=load("data/Vs22.0-beta-0.1-D_r-1.0-kappa-"+string(kappa)+"-rep-0.0-rate.txt");
J(2,:)=load("data/Vs22.0-beta-0.1-D_r-1.0-kappa-"+string(kappa)+"-rep-1.0-rate.txt");
J(3,:)=load("data/Vs22.0-beta-0.1-D_r-1.0-kappa-"+string(kappa)+"-rep-2.0-rate.txt");
J(4,:)=load("data/Vs22.0-beta-0.1-D_r-1.0-kappa-"+string(kappa)+"-rep-3.0-rate.txt");

errorbar(gam,mean(J,1),std(J,1),'x','Color',colour(i,:),LineWidth=2,MarkerSize=20,DisplayName='$\mathcal{V}_s='+string(Vs)+'\mu$m/s, $D_r$='+D_r_str...
    +'s$^{-1}$, $\beta$='+string(beta)+', $\kappa=$'+string(kappa));hold on
i=i+1
end


set(gca,'Xscale','log') 
set(gca,'Yscale','log')
xlabel('Fluid shear rate $\dot{\gamma} $ (s$^{-1}$)')
ylabel({"Net bacterial adhesion rate $J_A$ (cells/s$^{-1}$)"})
legend(location='southwest')
% Add scaling lines on top
loglog(gam(1:6),0.8e8*gam(1:6).^(1/3)./(2500*750),'k-.',HandleVisibility='off')
    loglog(gam(12:16),25e7./(2500*750)*gam(12:16).^(-1) ...
    ,'k-.',HandleVisibility='off')

title('(a)')

ylim([3,1e2])
  xlim([2e-2,6e1])
 exportgraphics(figure1,'flux_kappa.pdf')


%% Plot Fig s4b.

figure1=figure('units','inch','position',[0,0,7,6 ...
    ]);


figure
    J=[];

colour=[1 0.411764705882353 0.16078431372549];

J(1,:)=load("Vs10.0-beta-0.0-D_r-0.5-rep-0.0-rate.txt");
J(2,:)=load("Vs10.0-beta-0.0-D_r-0.5-rep-1.0-rate.txt");
J(3,:)=load("Vs10.0-beta-0.0-D_r-0.5-rep-2.0-rate.txt");
J(4,:)=load("Vs10.0-beta-0.0-D_r-0.5-rep-3.0-rate.txt");

errorbar(gam,mean(J,1),std(J,1),'.',Color=colour,LineWidth=2,MarkerSize=20,DisplayName='$\mathcal{V}_s='+string(Vs)+'\mu$m/s, $D_r$='+D_r_str...
    +'s$^{-1}$, $\beta$='+string(beta)+', $\alpha$=0');hold on

    J=[];
J(1,:)=load("Vs10.0-beta-0.0-D_r-0.5-rep-0.0-rate-dip.txt");
J(2,:)=load("Vs10.0-beta-0.0-D_r-0.5-rep-1.0-rate-dip.txt");
J(3,:)=load("Vs10.0-beta-0.0-D_r-0.5-rep-2.0-rate-dip.txt");
J(4,:)=load("Vs10.0-beta-0.0-D_r-0.5-rep-3.0-rate-dip.txt");



errorbar(gam,mean(J,1),std(J,1),'x',Color=colour,LineWidth=2,MarkerSize=20,DisplayName='$\mathcal{V}_s='+string(Vs)+'\mu$m/s, $D_r$='+D_r_str...
    +'s$^{-1}$, $\beta$='+string(beta)+', $\alpha$=31.82$\mu$m$^2$/s');hold on


    J=[];
J(1,:)=load("Vs22.0-beta-0.1-D_r-1.0-rep-0.0-rate.txt");
J(2,:)=load("Vs22.0-beta-0.1-D_r-1.0-rep-1.0-rate.txt");
J(3,:)=load("Vs22.0-beta-0.1-D_r-1.0-rep-2.0-rate.txt");
J(4,:)=load("Vs22.0-beta-0.1-D_r-1.0-rep-3.0-rate.txt");

errorbar(gam,mean(J,1),std(J,1),'r.',LineWidth=2,MarkerSize=20,DisplayName='$\mathcal{V}_s='+string(Vs)+'\mu$m/s, $D_r$='+D_r_str...
    +'s$^{-1}$, $\beta$='+string(beta)+', $\alpha$=0');hold on

    J=[];
J(1,:)=load("Vs22.0-beta-0.1-D_r-1.0-rep-0.0-rate-dip.txt");
J(2,:)=load("Vs22.0-beta-0.1-D_r-1.0-rep-1.0-rate-dip.txt");
J(3,:)=load("Vs22.0-beta-0.1-D_r-1.0-rep-2.0-rate-dip.txt");
J(4,:)=load("Vs22.0-beta-0.1-D_r-1.0-rep-3.0-rate-dip.txt");

errorbar(gam,mean(J,1),std(J,1),'rx',LineWidth=2,MarkerSize=20,DisplayName='$\mathcal{V}_s='+string(Vs)+'\mu$m/s, $D_r$='+D_r_str...
    +'s$^{-1}$, $\beta$='+string(beta)+', $\alpha$=31.82$\mu$m$^2$/s');hold on


    J=[];

        colour=[0.392156862745098 0.831372549019608 0.0745098039215686];

J(1,:)=load("Vs70.0-beta-0.0-D_r-2.0-rep-0.0-rate.txt");
J(2,:)=load("Vs70.0-beta-0.0-D_r-2.0-rep-1.0-rate.txt");
J(3,:)=load("Vs70.0-beta-0.0-D_r-2.0-rep-2.0-rate.txt");
J(4,:)=load("Vs70.0-beta-0.0-D_r-2.0-rep-3.0-rate.txt");

errorbar(gam,mean(J,1),std(J,1),'r.',Color=colour,LineWidth=2,MarkerSize=20,DisplayName='$\mathcal{V}_s='+string(Vs)+'\mu$m/s, $D_r$='+D_r_str...
    +'s$^{-1}$, $\beta$='+string(beta)+', $\alpha$=0');hold on

    J=[];
J(1,:)=load("Vs70.0-beta-0.0-D_r-2.0-rep-0.0-rate-dip.txt");
J(2,:)=load("Vs70.0-beta-0.0-D_r-2.0-rep-1.0-rate-dip.txt");
J(3,:)=load("Vs70.0-beta-0.0-D_r-2.0-rep-2.0-rate-dip.txt");
J(4,:)=load("Vs70.0-beta-0.0-D_r-2.0-rep-3.0-rate-dip.txt");

errorbar(gam,mean(J,1),std(J,1),'rx',Color=colour,LineWidth=2,MarkerSize=20,DisplayName='$\mathcal{V}_s='+string(Vs)+'\mu$m/s, $D_r$='+D_r_str...
    +'s$^{-1}$, $\beta$='+string(beta)+', $\alpha$=31.82$\mu$m$^2$/s');hold on


set(gca,'Xscale','log') 
set(gca,'Yscale','log')

xlabel('Fluid shear rate $\dot{\gamma} $ (s$^{-1}$)')
ylabel('Net bacterial adhesion rate $J_A$(cells/s$^{-1}$)')
legend(location='southwest')

ylim([3e-2,1e3])
  xlim([2e-2,6e1])
title('(b)')
exportgraphics(figure1,'flux-dipoles.pdf')
