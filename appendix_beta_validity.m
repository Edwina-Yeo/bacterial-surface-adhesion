%% Solving for the angular distribution using series solution from 
% J. Talbot, C. Antoine, Philippe Claudin, E. Somfai, T. Börzsönyi. Exploring noisy Jeffery orbits: A combined Fokker-Planck and Langevin analysis in two and three dimensions. Physical Review E , 2024, 110 (4), pp.044143. 10.1103/PhysRevE.110.044143 . hal-04794688
% close all

set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultTextFontname', 'CMU Serif')
set(groot,'DefaultAxesFontName', 'CMU Serif')
set(groot,'DefaultLineLineWidth',1.2)
set(groot,'DefaultAxesFontSize',20);
set(0, 'DefaultFigureRenderer', 'opengl');
jeff=@(beta,x)1./(-1/2 + beta/2*cos(2*x));
theta_vec=linspace(0,2*pi,100);


figure1=figure('units','inch','position',[0,0,15,6]);

beta_vec=linspace(0,1,30)
q_sol12=beta_vec
q_sol11=beta_vec
q_diff_exact11=beta_vec

q_diff_exact12=beta_vec
Q11_hinch=beta_vec
Q12_hinch=beta_vec
Q11_doi=beta_vec
Q12_doi=beta_vec
nexttile
per=0.1
for i=1:length(beta_vec)
Per=per;
beta=beta_vec(i);
sol=calc_expansion_sol(Per,beta,theta_vec);

q_sol12(i)= (beta .*Per)./((16+Per.^2)) ;
q_sol11(i) =  (beta .*Per^2)./(4 *(16+Per.^2));


Qyy=trapz(theta_vec,(sol).*(sin(theta_vec)).^2);
Qxx=trapz(theta_vec,(sol).*(cos(theta_vec)).^2);
Qxy=trapz(theta_vec,sol.*sin(theta_vec).*cos(theta_vec));

q_diff_exact11(i)=Qxx-1/2
q_diff_exact12(i)=Qxy


% Define: doi closure solution:
A = 1024*beta^3*Per^3 + 576*beta^3*Per^5 + 288*beta^5*Per^5;
B = -64*beta^2*Per^2 + 12*beta^2*Per^4 - 12*beta^4*Per^4;
C = sqrt(4*B^3 + A^2);
D = A + C;
D13 = D^(1/3);  


expr = (2*Per)/(3*beta) + (beta*Per)/3 ...
    - (1024*2^(2/3)*beta*Per)/(9*D^(2/3)) ...
    + (128*2^(2/3)*beta*Per^3)/(3*D^(2/3)) ...
    - (128*2^(2/3)*beta^3*Per^3)/(3*D^(2/3)) ...
    - (4*2^(2/3)*beta*Per^5)/(D^(2/3)) ...
    + (8*2^(2/3)*beta^3*Per^5)/(D^(2/3)) ...
    - (4*2^(2/3)*beta^5*Per^5)/(D^(2/3)) ...
    + (128*2^(1/3))/(9*D^(1/3)) ...
    - (8*2^(1/3)*Per^2)/(3*D^(1/3)) ...
    + (8*2^(1/3)*beta^2*Per^2)/(3*D^(1/3)) ...
    + (2^(2/3)*D^(1/3))/(9*beta^2*Per^2) ...
    - D^(2/3)/(36*2^(2/3)*beta^3*Per^3);

Q11_doi(i) = expr / (2*Per);
Q12_doi(i) = -4/(3*beta*Per) ...
    - B/(6*2^(2/3)*beta^2*Per^2*D13) ...
    + D13/(12*2^(1/3)*beta^2*Per^2);


b=beta;
inner_sqrt = sqrt(64 * b^2 * Per^2 + (16 + Per^2 - b^2 * Per^2)^2);
outer_sqrt = sqrt(16 - Per^2 + b^2 * Per^2 + inner_sqrt);



Q12_hinch(i) = (-4 + outer_sqrt / sqrt(2)) / (2 * b * Per);

term1 = (2 * Per) / b;
    
inner_sqrt = sqrt(64 * b^2 * Per^2 + (16 + Per^2 - b^2 * Per^2)^2);
outer_expr = 16 - Per^2 + b^2 * Per^2 + inner_sqrt;
outer_sqrt = sqrt(outer_expr);

term2 = (4 * sqrt(2) * outer_sqrt) / (b * Per);
term3 = (Per * outer_sqrt) / (2 * sqrt(2) * b);
term4 = (b * Per * outer_sqrt) / (2 * sqrt(2));
term5 = (outer_expr)^(3/2) / (4 * sqrt(2) * b * Per);

Q11_hinch(i) = (term1 + term2 - term3 + term4 - term5) / (4 * Per);

end

subplot(1,2,1)
plot(beta_vec,q_diff_exact11,'k.',MarkerSize=20,DisplayName='(S33), $Pe_r=0.1$'); hold on 
plot(beta_vec,q_sol11,'k-',DisplayName='Continuum model, $Pe_r=0.1$'); hold on 
% plot(beta_vec,Q11_hinch,'k--',DisplayName='Hinch \& Leal Closure, $Pe_r=0.1$'); hold on
% plot(beta_vec,Q11_doi,'k-.',DisplayName='Doi et al. closure, $Pe_r=0.1$'); hold on
subplot(1,2,2)
plot(beta_vec,q_diff_exact12,'k.',MarkerSize=20,DisplayName='(S33), $Pe_r=0.1$'); hold on 
plot(beta_vec,q_sol12,'k-',DisplayName='Continuum model, $Pe_r=0.1$'); hold on
% plot(beta_vec,Q12_hinch,'k--',DisplayName='Hinch \& Leal Closure, $Pe_r=0.1$'); hold on
% plot(beta_vec,Q12_doi,'k-',DisplayName='Doi et al. closure, $Pe_r=0.1$'); hold on

per=5
for i=1:length(beta_vec)
Per=per;
beta=beta_vec(i);
sol=calc_expansion_sol(Per,beta,theta_vec);

q_sol12(i)= (beta .*Per)./((16+Per.^2)) ;
q_sol11(i) =  (beta .*Per^2)./(4 *(16+Per.^2));


Qyy=trapz(theta_vec,(sol).*(sin(theta_vec)).^2);
Qxx=trapz(theta_vec,(sol).*(cos(theta_vec)).^2);
Qxy=trapz(theta_vec,sol.*sin(theta_vec).*cos(theta_vec));

q_diff_exact11(i)=Qxx-1/2
q_diff_exact12(i)=Qxy


% Define: doi closure solution:
A = 1024*beta^3*Per^3 + 576*beta^3*Per^5 + 288*beta^5*Per^5;
B = -64*beta^2*Per^2 + 12*beta^2*Per^4 - 12*beta^4*Per^4;
C = sqrt(4*B^3 + A^2);
D = A + C;
D13 = D^(1/3);  


expr = (2*Per)/(3*beta) + (beta*Per)/3 ...
    - (1024*2^(2/3)*beta*Per)/(9*D^(2/3)) ...
    + (128*2^(2/3)*beta*Per^3)/(3*D^(2/3)) ...
    - (128*2^(2/3)*beta^3*Per^3)/(3*D^(2/3)) ...
    - (4*2^(2/3)*beta*Per^5)/(D^(2/3)) ...
    + (8*2^(2/3)*beta^3*Per^5)/(D^(2/3)) ...
    - (4*2^(2/3)*beta^5*Per^5)/(D^(2/3)) ...
    + (128*2^(1/3))/(9*D^(1/3)) ...
    - (8*2^(1/3)*Per^2)/(3*D^(1/3)) ...
    + (8*2^(1/3)*beta^2*Per^2)/(3*D^(1/3)) ...
    + (2^(2/3)*D^(1/3))/(9*beta^2*Per^2) ...
    - D^(2/3)/(36*2^(2/3)*beta^3*Per^3);

Q11_doi(i) = expr / (2*Per);
Q12_doi(i) = -4/(3*beta*Per) ...
    - B/(6*2^(2/3)*beta^2*Per^2*D13) ...
    + D13/(12*2^(1/3)*beta^2*Per^2);


b=beta;
inner_sqrt = sqrt(64 * b^2 * Per^2 + (16 + Per^2 - b^2 * Per^2)^2);
outer_sqrt = sqrt(16 - Per^2 + b^2 * Per^2 + inner_sqrt);



Q12_hinch(i) = (-4 + outer_sqrt / sqrt(2)) / (2 * b * Per);

term1 = (2 * Per) / b;
    
inner_sqrt = sqrt(64 * b^2 * Per^2 + (16 + Per^2 - b^2 * Per^2)^2);
outer_expr = 16 - Per^2 + b^2 * Per^2 + inner_sqrt;
outer_sqrt = sqrt(outer_expr);

term2 = (4 * sqrt(2) * outer_sqrt) / (b * Per);
term3 = (Per * outer_sqrt) / (2 * sqrt(2) * b);
term4 = (b * Per * outer_sqrt) / (2 * sqrt(2));
term5 = (outer_expr)^(3/2) / (4 * sqrt(2) * b * Per);

Q11_hinch(i) = (term1 + term2 - term3 + term4 - term5) / (4 * Per);

end

subplot(1,2,1)
plot(beta_vec,q_diff_exact11,'b.',MarkerSize=20,DisplayName='(S33), $Pe_r=5$'); hold on 
plot(beta_vec,q_sol11,'b-',DisplayName='Continuum model, $Pe_r=5$'); hold on 
% plot(beta_vec,Q11_hinch,'b--',DisplayName='Hinch \& Leal Closure, $Pe_r=0.1$'); hold on
% plot(beta_vec,Q11_doi,'b-.',DisplayName='Doi et al. closure, $Pe_r=0.1$'); hold on
subplot(1,2,2)
plot(beta_vec,q_diff_exact12,'b.',MarkerSize=20,DisplayName='(S33), $Pe_r=5$'); hold on 
plot(beta_vec,q_sol12,'b-',DisplayName='Continuum model, $Pe_r=5$'); hold on 
% plot(beta_vec,Q12_hinch,'b--',DisplayName='Hinch \& Leal Closure, $Pe_r=0.1$'); hold on
% plot(beta_vec,Q12_doi,'b-.',DisplayName='Doi et al. closure, $Pe_r=0.1$'); hold on

per=100
for i=1:length(beta_vec)
Per=per;
beta=beta_vec(i);
sol=calc_expansion_sol(Per,beta,theta_vec);

q_sol12(i)= (beta .*Per)./((16+Per.^2)) ;
q_sol11(i) =  (beta .*Per^2)./(4 *(16+Per.^2));


Qyy=trapz(theta_vec,(sol).*(sin(theta_vec)).^2);
Qxx=trapz(theta_vec,(sol).*(cos(theta_vec)).^2);
Qxy=trapz(theta_vec,sol.*sin(theta_vec).*cos(theta_vec));

q_diff_exact11(i)=Qxx-1/2
q_diff_exact12(i)=Qxy

% Define: doi closure solution:
A = 1024*beta^3*Per^3 + 576*beta^3*Per^5 + 288*beta^5*Per^5;
B = -64*beta^2*Per^2 + 12*beta^2*Per^4 - 12*beta^4*Per^4;
C = sqrt(4*B^3 + A^2);
D = A + C;
D13 = D^(1/3);  


expr = (2*Per)/(3*beta) + (beta*Per)/3 ...
    - (1024*2^(2/3)*beta*Per)/(9*D^(2/3)) ...
    + (128*2^(2/3)*beta*Per^3)/(3*D^(2/3)) ...
    - (128*2^(2/3)*beta^3*Per^3)/(3*D^(2/3)) ...
    - (4*2^(2/3)*beta*Per^5)/(D^(2/3)) ...
    + (8*2^(2/3)*beta^3*Per^5)/(D^(2/3)) ...
    - (4*2^(2/3)*beta^5*Per^5)/(D^(2/3)) ...
    + (128*2^(1/3))/(9*D^(1/3)) ...
    - (8*2^(1/3)*Per^2)/(3*D^(1/3)) ...
    + (8*2^(1/3)*beta^2*Per^2)/(3*D^(1/3)) ...
    + (2^(2/3)*D^(1/3))/(9*beta^2*Per^2) ...
    - D^(2/3)/(36*2^(2/3)*beta^3*Per^3);

Q11_doi(i) = expr / (2*Per);
Q12_doi(i) = -4/(3*beta*Per) ...
    - B/(6*2^(2/3)*beta^2*Per^2*D13) ...
    + D13/(12*2^(1/3)*beta^2*Per^2);


b=beta;
inner_sqrt = sqrt(64 * b^2 * Per^2 + (16 + Per^2 - b^2 * Per^2)^2);
outer_sqrt = sqrt(16 - Per^2 + b^2 * Per^2 + inner_sqrt);



Q12_hinch(i) = (-4 + outer_sqrt / sqrt(2)) / (2 * b * Per);

term1 = (2 * Per) / b;
    
inner_sqrt = sqrt(64 * b^2 * Per^2 + (16 + Per^2 - b^2 * Per^2)^2);
outer_expr = 16 - Per^2 + b^2 * Per^2 + inner_sqrt;
outer_sqrt = sqrt(outer_expr);

term2 = (4 * sqrt(2) * outer_sqrt) / (b * Per);
term3 = (Per * outer_sqrt) / (2 * sqrt(2) * b);
term4 = (b * Per * outer_sqrt) / (2 * sqrt(2));
term5 = (outer_expr)^(3/2) / (4 * sqrt(2) * b * Per);

Q11_hinch(i) = (term1 + term2 - term3 + term4 - term5) / (4 * Per);

end

subplot(1,2,1)
plot(beta_vec,q_diff_exact11,'r.',MarkerSize=20,DisplayName='(S33), $Pe_r=100$'); hold on 
plot(beta_vec,q_sol11,'r-',DisplayName='(S7-8) closure, $Pe_r=100$'); hold on 
% plot(beta_vec,Q11_hinch,'r--',DisplayName='HL1 closure, $Pe_r=0.1$'); hold on
% plot(beta_vec,Q11_doi,'r-.',DisplayName='Doi closure, $Pe_r=0.1$'); hold on
legend(Location='northwest')
title('(a)')
xlabel('Bretherton parameter, $\beta$')
ylabel('$Q_{xx}$')
ylim([0,0.45])
subplot(1,2,2)
plot(beta_vec,q_diff_exact12,'r.',MarkerSize=20,DisplayName='(S33), $Pe_r=5$'); hold on 
plot(beta_vec,q_sol12,'r-',DisplayName='Continuum model, $Pe_r=5$'); hold on
% plot(beta_vec,Q12_hinch,'r--',DisplayName='HL1 closure, $Pe_r=0.1$'); hold on
% plot(beta_vec,Q12_doi,'r-.',DisplayName='Doi closure, $Pe_r=0.1$'); hold on
title('(b)')
xlabel('Bretherton parameter, $\beta$')
ylabel('$Q_{xy}$')
ylim([0,0.25])
exportgraphics(figure1,'beta_validity.pdf')
%



function sol=calc_expansion_sol(Per,beta,theta_vec)
N=30;% number of terms from the expansion
% beta=0.5;
% Per=1;
sol_vec=zeros(2*N);% (a=1:N),(b=1+N:2N)
a0=1/(2*pi);
% Create matrix for solving

A=zeros(2*N,2*N);
% create a rows first 

B=zeros(2*N,1);
% -2a1+beta a_{2})-8/Per*bk==- 2beta a0

A(1,1)=-2;
B(1)=-2*beta*a0; % move the a0 term to RHS
A(1,2)=beta;
A(1,N+1)=-8/Per;

% -2ak+beta(a{k-1}+a_{k+1})-8/Per*k*bk==0


for i=2:N-1
A(i,i)=-2;
A(i,i-1)=beta;
A(i,i+1)=beta;
A(i,i+N)=-8/Per*i;
end
% -2aN+beta(aN-1)-8/Per*k*bk==0

A(N,N)=-2;
A(N,N-1)=beta;
A(N,2*N)=-8/Per*N;


% 2b1-beta(b_2})-8/Per*a1==0
A(N+1,N+1)=2;
A(N+1,N+2)=-beta;
A(N+1,1)=-8/Per;
% 2bk-beta(b{k-1}+b_{k+1})-8/Per*k*ak==0
for i=N+2:2*N-1
A(i,i)=2;
A(i,i-1)=-beta;
A(i,i+1)=-beta;
A(i,i-N)=-8/Per*(i-N);

end

% 2bN-beta(b{N-1})-8/Per*k*aN==0
A(2*N,2*N)=2;
A(2*N,2*N-1)=-beta;
A(2*N,N)=-8/Per*(N);


coef_vec=A\B;
% theta_vec=linspace(0,2*pi,100);
sol=a0*ones(size(theta_vec));
for i=1:N
sol=sol+coef_vec(i)*cos(2*i*theta_vec);
end

for i=N+1:2*N
sol=sol+coef_vec(i)*sin(2*(i-N)*theta_vec);

end


end

