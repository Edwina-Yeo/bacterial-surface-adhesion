
% Generates Fig S

Pe_eff_2D=@(beta,Vss,Pers)  (16 + (4 - beta^2)*Pers.^2).*Vss.^2 ./(32 + (2- beta)*(1- beta)*Pers.^2).*(16+Pers.^2)./(4*Vss.^2.*Pers);
% Pe_eff_3D=@(b,Vss,per) -(4*per.*(-44100 + (-1225 + 3*(980 - 451*b)*b)*per.^2))./(3*(-1764 + (-49 + 3*b^2)*per.^2).*(-100 + (-25 + 9*b.^2)*per.^2));
  
Pe_eff_3D=@(b,Vss,per) (4*per.*(88200 + (2450 + 3*b*(-1225 + 461*b)).*per.^2))./(3*(-1764 + (-49 + 3*b^2).*per.^2).*(-400 + (-25 + 9*b^2).*per.^2))
close all

    x_int=3.12013; % integral of x^(-1/3) between 0.5 and 3
coeff_line=x_int*3^(1/3)/gamma(1/3);

pers=logspace(-2,2,100)
figure1=figure('units','inch','position',[0,0,8 ...
    ,6 ...
    ]);
% subplot(1,2,1)
semilogx(pers,coeff_line*Pe_eff_2D(0,1,pers).^(-1/3),'k'); hold on
semilogx(pers,coeff_line*Pe_eff_2D(0.5,1,pers).^(-1/3),'k-.'); hold on


semilogx(pers,coeff_line*Pe_eff_3D(0,1,pers).^(1/3),'b'); hold on
semilogx(pers,coeff_line*Pe_eff_3D(0.5,1,pers).^(1/3),'b-.'); hold on

xlabel("Rotational P\'eclet number $Pe_r$")
ylabel('Net bacterial adhesion rate, $\bar{J}$')
lk=legend('2D - spherical, $\beta=0$', '2D - elongated, $\beta=0.5$','3D - spherical, $\beta=0$', '3D - elongated, $\beta=0.5$')

lk.Location='south'
% title('(a)')

% subplot(1,2,2)
% 
% semilogx(pers,Pe_eff_2D(0.5,1,pers).^(-1/3)./Pe_eff_2D(0,1,pers).^(-1/3),'k-.'); hold on
% semilogx(pers,Pe_eff_3D(0.5,1,pers).^(1/3)./Pe_eff_3D(0,1,pers).^(1/3) ...
%     ,'b-.'); hold on
% 
% xlabel("Rotational P\'eclet number $Pe_r$")
% ylabel('Relative adhesion rate')
% lk=legend('2D, $\beta=0.5$', '3D,  $\beta=0.5$')
% title('(b)')

exportgraphics(figure1,'3D.pdf')

% semilogx(pers,((4+pers.^2)./(2*pers)).^(-1/3)); hold on
