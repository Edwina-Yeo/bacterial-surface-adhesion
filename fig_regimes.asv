% Creates Fig 2d. 

set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultTextFontname', 'CMU Serif')
set(groot,'DefaultAxesFontName', 'CMU Serif')
set(groot,'DefaultLineLineWidth',3)
set(groot,'DefaultAxesFontSize',30);
set(0, 'DefaultFigureRenderer', 'opengl');

vs=logspace(-10,0,1000);%(eps^4,1,1000);
per=logspace(-7,5,1000);
[V,P]=meshgrid(vs,per);

% Using the spherical diffusivity. 
eps=(2*V.^2.*P./((4+P.^2))).^(1/3); % Epsilon over the 2D domain
Peff=(2*V.^2.*P./((4+P.^2)));% Epsilon over the 2D domain

% Find edges of domains from contour plots
% Region B:
body=5; %microns
L=10000; %microns  1cm 
delta=body/L;
kk1=contour(V,P,log(eps)-log(delta),[0,1.1]);

% Region C:
kk2=contour(V,P,log(1./eps)-log(P),[0,1.1]);

% Region D:
kk3=contour(V,P,log(eps)-log(V),[0,1.1]);


close all
figure1=figure('units','inch','position',[0,0,11,9]);

% Region B:
p=patch([kk1(1,2:10:1570),xmin,xmin,kk1(1,2)],[kk1(2,2:10:1570),ymin,ymax,kk1(2,2)],[0.6,0.6,0.7]); hold on
p.HandleVisibility='off';
p.EdgeColor='none';
plot(kk1(1,222:1570),kk1(2,222:1570),'k-.',LineWidth=2,HandleVisibility='off');hold on

% Region D:
plot(kk3(1,10:1300),kk3(2,10:1300),'k-.',LineWidth=2,HandleVisibility='off');hold on
p=patch([kk3(1,100:10:1500),1,1,kk3(1,100)],[kk3(2,100:10:1500),ymax,ymin,kk3(2,100)],[0.7,0.7,0.6]); hold on
p.HandleVisibility='off';
p.EdgeColor='none';

% Region C:
p=patch([xmax,kk2(1,2:10:800),1,1],[ymax,kk2(2,2:10:800 ),ymax,kk2(2,10)],[0.7,0.7,0.7]); hold on
p.HandleVisibility='off';
p.EdgeColor='none';
plot(kk2(1,2:10:800),kk2(2,2:10:800),'k-.',LineWidth=2,HandleVisibility='off')

% Region A:
plot([1,1],[ymin,ymax],'k-.',LineWidth=2,HandleVisibility='off');hold on
p=patch([1 100 100 1], [1e-6 1e-6 ymax ymax], [0.8,0.8,0.8]); hold on
p.HandleVisibility='off';
p.EdgeColor='none';

% Plot the effective Peclet number 
k1=surf(V,P,log(Peff),HandleVisibility='off'); hold on
cmap=load('CustomCmap.mat').CustomColormap;
colormap(cmap)
k1.EdgeColor='none';
colorbar()
view(2)
clim([-22,-2])

% Relabel xticks as we are plotting Log(Pe_eff)
colorbar('TickLabelInterpreter','latex',...
    'Direction','reverse','Ticks',[-40 -35 -30 -25 -20 -15 -10 -5 -2],...
    'TickLabels',{'$10^{40}$','$10^{35}$','$10^{30}$','$10^{25}$','$10^{20}$','$10^{15}$','$10^{10}$','$10^{5}$','0.01'},'TickLabelInterpreter','latex',...
    'Direction','reverse');

xlabel('Relative swimming speed, $V_s$')
ylabel('Rotational Peclet number $Pe_r$')
legend(Location="southwest")

% axes limits. 
xmin=5e-6;
xmax=10;
ymin=5e-4;
ymax=1e4;

ylim([ymin,ymax])
xlim([xmin,xmax])

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')


% Add data from bacteria in applications-----------------------

% Arteries
L=2600;%2.6mm
gamma_min_art=100 ;
gamma_max_art=350;
gamma_art=100/2+350/2;

[vs,vs_min,vs_max,per,per_min,per_max]=    ecoli_params(gamma_art,gamma_min_art,gamma_max_art,L);
errorbar(vs,per,per-per_min ...
    ,per_max-per,vs-vs_min ...
    ,-vs+vs_max,'c.',LineWidth=2,MarkerSize=50,DisplayName='Artery') 
[vs,vs_min,vs_max,per,per_min,per_max]=    PA_params(gamma_art,gamma_min_art,gamma_max_art,L);

errorbar(vs,per,per-per_min ...
    ,per_max-per,vs-vs_min ...
    ,-vs+vs_max,'cs',LineWidth=2,MarkerSize=18,MarkerFaceColor='c',DisplayName='Artery, {\it P. Aeruginosa}',HandleVisibility='off') 

% River bed flow
L=10000; %microns  1cm 
gamma_min_art=0.0059/0.001; % 1/s: 0,1Lmin /2mm catheter 
gamma_max_art=0.06/0.001;% 1/s: 4Lmin /0.5mm catheter 
gamma_art=0.0217/0.001;

col=[0.07,0.62,1]; 
[vs,vs_min,vs_max,per,per_min,per_max]=    ecoli_params(gamma_art, ...
    gamma_min_art,gamma_max_art,L);
errorbar(vs,per,per-per_min ...
    ,per_max-per,vs-vs_min ...
    ,-vs+vs_max,'.',MarkerEdgeColor=col,Color=col,LineWidth=2,MarkerSize=50,DisplayName='River bed') 
[vs,vs_min,vs_max,per,per_min,per_max]=    PA_params(gamma_art, ...
    gamma_min_art,gamma_max_art,L);
errorbar(vs,per,per-per_min ...
    ,per_max-per,vs-vs_min ...
    ,-vs+vs_max,'s',Color=col,LineWidth=2,MarkerSize=18,MarkerFaceColor=col,DisplayName='River bed, {\it P. Aeruginosa}',HandleVisibility='off') 

% Human Gut 
gamma_min_gut=0.0002/1e-3;
gamma_max_gut=0.008/1e-3;
gamma_gut=gamma_min_gut/2+gamma_max_gut/2;
L=12500; %microns  1.25cm 


[vs,vs_min,vs_max,per,per_min,per_max]=    ecoli_params(gamma_gut,gamma_min_gut,gamma_max_gut,L);
errorbar(vs,per,per-per_min ...
    ,per_max-per,vs-vs_min ...
    ,-vs+vs_max,'g.',LineWidth=2,MarkerSize=50,DisplayName='Gut'); 

[vs,vs_min,vs_max,per,per_min,per_max]=    PA_params(gamma_gut,gamma_min_gut,gamma_max_gut,L);
errorbar(vs,per,per-per_min ...
    ,per_max-per,vs-vs_min ...
    ,-vs+vs_max,'gs',LineWidth=2,MarkerSize=18,MarkerFaceColor='g',DisplayName='Gut, {\it P. Aeruginosa}',HandleVisibility='off'); 


% Small catheter
L=1000; %microns  (1mm) 

gamma_min=4*0.1/(pi*1^3); % 1/s: 0,1Lmin /1mm catheter 
gamma_max_cath=4*4/(pi*1^3);% 1/s: 4Lmin /1mm catheter 
gamma_cath=gamma_max_cath/2+gamma_min/2;

[vs,vs_min,vs_max,per,per_min,per_max]=    ecoli_params(gamma_cath,gamma_min,gamma_max_cath,L);
errorbar(vs,per,per-per_min ...
    ,per_max-per,vs-vs_min ...
    ,-vs+vs_max,'r.',LineWidth=2,MarkerSize=50,DisplayName='Small catheter') ;
[vs,vs_min,vs_max,per,per_min,per_max]=    PA_params(gamma_cath,gamma_min,gamma_max_cath,L);
errorbar(vs,per,per-per_min ...
    ,per_max-per,vs-vs_min ...
    ,-vs+vs_max,'rs',LineWidth=2,MarkerSize=18,MarkerFaceColor='r',DisplayName='Small catheter, {\it P. Aeruginosa}',HandleVisibility='off') 


% Large catheter
L=1500; % 1.5mm
gamma_min=4*0.1/(pi*1.5^3); % 1/s: 0,1Lmin /1.5mm catheter 
gamma_max_cath=4*4/(pi*1.5^3);% 1/s: 4Lmin /1.5m catheter 
gamma_cath=gamma_max_cath/2+gamma_min/2;


[vs,vs_min,vs_max,per,per_min,per_max]=    ecoli_params(gamma_cath,gamma_min,gamma_max_cath,L);
errorbar(vs,per,per-per_min ...
    ,per_max-per,vs-vs_min ...
    ,-vs+vs_max,'k.',LineWidth=2,MarkerSize=50,DisplayName='Large catheter')
[vs,vs_min,vs_max,per,per_min,per_max]=    PA_params(gamma_cath,gamma_min,gamma_max_cath,L);
errorbar(vs,per,per-per_min ...
    ,per_max-per,vs-vs_min ...
    ,-vs+vs_max,'ks',LineWidth=2,MarkerSize=18,MarkerFaceColor='k',DisplayName='Large catheter, {\it P. Aeruginosa}',HandleVisibility='off') 

% Add region labels 
annotation(figure1,'textbox',...
    [0.570917508417511 0.806972001763665 0.040316458101623 0.0912698412698414],...
    'String','C: Horizontal Swimming',...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',27,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(figure1,'textbox',...
    [0.757952441077441 0.792660010240652 0.0802881150018602 0.105990783410138],...
    'String','A: Weak flow',...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',27,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(figure1,'textbox',...
    [0.150936447811448 0.78832138590203 0.0802881150018606 0.105990783410138],...
    'String','B: Passive tracers',...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',27,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(figure1,'textbox',...
    [0.568865740740744 0.18616089207487 0.0802881150018597 0.105990783410138],...
    'String',{'D:','Flow-bacterial','coupling'},...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',27,...
    'FitBoxToText','off',...
    'EdgeColor','none');

exportgraphics(figure1,'reigimes-combined.jpg',"Resolution",600)



% takes shear rate in (1/s) and length in microns and returns model params
% for E coli
function [vs,vs_min,vs_max,per,per_min,per_max]=ecoli_params(gamma,gamma_min, gamma_max,L)

% L_u is the velocity scale thats the thickenss of the domain in the y direction:NB not boundary layer 

U=gamma*L;
U_min=gamma_min*L;
U_max=gamma_max*L;
L_model=L;

dr=1; % estimated using tumbling rate 

V=22.00;% micron/s
V_min=22-5;
V_max=22+5;

vs=V/U;
vs_min=V_min/U_max;
vs_max=V_max/U_min;
per=U/(dr*L_model);
per_min=U_min/(dr*L_model);
per_max=U_max/(dr*L_model);


end
% for P.aeruginosa
function [vs,vs_min,vs_max,per,per_min,per_max]=PA_params(gamma,gamma_min, gamma_max,L)

% L_u is the velocity scale thats the thickenss of the domain in the y direction:NB not boundary layer 

U=gamma*L;
U_min=gamma_min*L;
U_max=gamma_max*L;
L_model=L;

dr=0.036; %(rad^2/s) no range given. 


V=22.00;% micron/s
V_min=22-6;
V_max=22+6;

vs=V/U;
vs_min=V_min/U_max;
vs_max=V_max/U_min;
per=U/(dr*L_model);
per_min=U_min/(dr*L_model);
per_max=U_max/(dr*L_model);

end
