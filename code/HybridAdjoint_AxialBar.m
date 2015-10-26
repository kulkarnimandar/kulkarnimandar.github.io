clear all; close all;
% Axial bar parameters
E = 110.3e9; % titanium in Pa
L = 5; % in meters
A = .2*.2; % in m^2
f_0 = 1; % N/m
% Number of elements
N = 50;

% All elements of equal lengths
Le = L/N;

% Elemental and global stiffness matrix
Ke = E*A/Le*[1 -1;-1 1];
K = zeros(N+1,N+1);
for i=1:N
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)+Ke;
end

% Elemental and global forces
F = zeros(N+1,1);
Fprime1 = zeros(N+1,1);
Fprime2 = zeros(N+1,1);
for i = 1:N
    xmid = (i-1/2)*Le;
    Fmid = f_0*(L-xmid);
    Fe = [Fmid*Le/2;Fmid*Le/2];
    F(i:i+1,1) = F(i:i+1,1) + Fe;
    %
    %Parameterization 1: fprime = 1
    F_prime_mid_param1 = f_0;
    Fe_prime_param1 = [F_prime_mid_param1*Le/2;F_prime_mid_param1*Le/2];
    Fprime1(i:i+1,1) = Fprime1(i:i+1,1) + Fe_prime_param1;
    %
    %Parameterization 2: fprime = 0
    F_prime_mid_param2 = 0;
    Fe_prime_param2 = [F_prime_mid_param2*Le/2;F_prime_mid_param2*Le/2];
    Fprime2(i:i+1,1) = Fprime2(i:i+1,1) + Fe_prime_param2;
end

% Partitioning stiffness matrix
cc = 1;
uu = setdiff(1:N+1,cc);
%
K_uu = K(uu,uu);
K_uc = K(uu,cc);

% Essential (or Geometric) BCs
u_cc = 0; % Geometric BC for u is u(0) = 0;
uprime1_cc = 0; % Geometric BC for uprime1 is uprime1(0) = 0;
uprime2_cc = L^2/(2*A*E); % Geometric BC for uprime2 is uprime2(0) = L^2/(2*A*E);

% Unconstrained Loading
F_uu = F(uu,1) - K_uc*u_cc;
Fprime1_uu = Fprime1(uu,1) - K_uc*uprime1_cc;
Fprime2_uu = Fprime2(uu,1) - K_uc*uprime2_cc;

% Solve for u
u_uu = K_uu\F_uu;
uprime1_uu = K_uu\Fprime1_uu;
uprime2_uu = K_uu\Fprime2_uu;

x = (linspace(0,L,N+1))';
u = [u_cc;u_uu]; 
uprime1 = [uprime1_cc;uprime1_uu]; 
uprime2 = [uprime2_cc;uprime2_uu]; 

% Spatial Gradients and design velocities
grad_u = (x.^2/2 - L*x + L^2/2)/(A*E);
V1 = x/L;
V2 = x/L - 1;

% Total sensitivities
conv1 = grad_u.*V1;
conv2 = grad_u.*V2;
udot1 = uprime1 + conv1;
udot2 = uprime2 + conv2;

% Non-dimensionalize
u_nd = u/(L^3/(6*E*A));
udot1_nd = udot1/(L^2/(2*A*E));
udot2_nd = udot2/(L^2/(2*A*E));
uprime1_nd = uprime1/(L^2/(2*A*E));
uprime2_nd = uprime2/(L^2/(2*A*E));
conv1_nd = conv1/(L^2/(2*A*E));
conv2_nd = conv2/(L^2/(2*A*E));

% Exact solutions
xex = linspace(0,L,101);
uex = f_0*(3*L^2*xex - 3*L*xex.^2 + xex.^3)/(6*E*A);
uexnd = uex/(f_0*L^3/(6*E*A));
%
udotex = (xex.^3/(2*L) - 3*xex.^2/2 + 3*L*xex/2)/(A*E);
udotex_nd = udotex/(L^2/(2*A*E));
%
uprimeex = (-xex.^2/2 + L*xex)/(A*E);
uprimeex_nd = uprimeex/(L^2/(2*A*E));

% Plot FE and exact solutions
set(0,'defaultlinelinewidth',2)
set(0,'defaultaxesfontsize',13)
set(0,'defaultaxesfontname','Times')

figure(1)
plot(xex,uexnd,'r');
hold on;
plot(x,u_nd,'bo');
xlabel('x [m]','interpreter','latex');
ylabel('$u/u_{tip}$ [m/m]','interpreter','latex');
% title('Displacement of axial bar with distributed load');
legend('Analytic solution (exact)','Finite Element solution','location','southeast')
axis([0 L 0 1.1])
grid on;

figure(2)
% plot(xex,uprimeex_nd,'r');
% hold on;
plot(x,uprime1_nd,'+');%,'markersize',12);
hold on;
plot(x,uprime2_nd,'r^');
% plot(x,[conv1_nd,conv2_nd],'*');
xlabel('x [m]','interpreter','latex');
ylabel('$u^{\prime}/\dot{u}_{tip}$','interpreter','latex');
% title('Local sensitivity of displacement to length of axial bar');

figure(3)
plot(xex,udotex_nd,'k');
hold on;
plot(x,udot1_nd,'+');%,'markersize',12);
plot(x,udot2_nd,'r^');%,'markersize',9);
xlabel('x [m]','interpreter','latex');
ylabel('$\dot{u}/\dot{u}_{tip}$','interpreter','latex');
% title('Total sensitivity of displacement to length of axial bar');

%% Calculation of local sensitivity by adjoint formulation
% g = z'*u
% g_prime = z'*u_prime + lambda'*(K_uu*u_prime - f_prime)
% Set (z' + lambda'*K_uu = 0), Thus solve for lambda: lambda = -(K_uu')^-1*z
% g_prime = -lambda'*f_prime
% 
% Here f_prime will be calculated for each design variable, but the
% inversion of K_uu matrix is done only once!

NodeMat = 2:floor(N/5):N+1;
gprime1 = zeros(length(NodeMat),1);
gprime2 = zeros(length(NodeMat),1);
lambdaMat = zeros(length(uu),length(NodeMat));

for j = 1:length(NodeMat)
    Node = NodeMat(j);
    z = zeros(size(u));
    z(Node) = 1;
    z_uu = z(uu);
    lambda = -K_uu'\z_uu;
    lambdaMat(:,j) = lambda;
    %
    gprime1(j) = - lambda'*Fprime1_uu;
    gprime2(j) = - lambda'*Fprime2_uu;
end
%
gdot1 = gprime1 + conv1(NodeMat);
gdot2 = gprime2 + conv2(NodeMat);
%
gprime1_nd = gprime1/(L^2/(2*A*E));
gprime2_nd = gprime2/(L^2/(2*A*E));
gdot1_nd = gdot1/(L^2/(2*A*E));
gdot2_nd = gdot2/(L^2/(2*A*E));

figure(2)
plot(x(NodeMat),gprime1_nd,'o','MarkerSize',12);
plot(x(NodeMat),gprime2_nd,'rsq','MarkerSize',12);
legend('Direct CSA, Param. 1','Direct CSA, Param. 2','Hybrid Adj., Param. 1','Hybrid Adj., Param. 2','location','southeast')
axis([0 L 0 1.1]);
grid on;

figure(3)
plot(x(NodeMat),gdot1_nd,'o','MarkerSize',12);
plot(x(NodeMat),gdot2_nd,'rsq','MarkerSize',12);
legend('Exact solution','Direct CSA, Param. 1','Direct CSA, Param. 2','Hybrid Adj., Param. 1','Hybrid Adj., Param. 2','location','southeast')
grid on;

figure()
plot(x,conv1_nd,'b.','MarkerSize',18);
hold on;
plot(x,conv2_nd,'rx','MarkerSize',10);
xlabel('x [m]','interpreter','latex');
ylabel('$[\nabla_xu.V] /\dot{u}_{tip}$','interpreter','latex');
legend('Convective term, Param. 1','Convective term, Param. 2','location','southeast')
grid on;

lambda_nd = lambdaMat/(1/(E*A/L));
figure()
% h(1)=plot(x(uu),Fprime1_uu,'bo');%,'MarkerSize',10);
% h(2)=plot(x(uu),Fprime1(uu),'b+','MarkerSize',10);
% plot(x(uu),Fprime2_uu,'ro','MarkerSize',10);
% plot(x(uu),Fprime2(uu),'rx','MarkerSize',10);
h(3)=plot(x(uu),lambda_nd(:,1),'g+');
hold on;
h(4)=plot(x(uu),lambda_nd(:,2),'cx');
h(5)=plot(x(uu),lambda_nd(:,3),'m*');
h(6)=plot(x(uu),lambda_nd(:,4),'ro');
h(7)=plot(x(uu),lambda_nd(:,5),'kd');
xlabel('x [m]','interpreter','latex');
ylabel('Loads and the adjoint variable');
grid on;
axis([-inf inf -1.2 0.2]);
% legend('F^{\prime}_{Local} Param. 1','F^{\prime}_{1} Param. 1','F^{\prime}_{Local} Param. 2','F^{\prime}_{1} Param. 2','\lambda','location','southeast')
% legend('F^{\prime}_{Local} Param. 1','F^{\prime}_{1} Param. 1','\lambda_1','\lambda_2','\lambda_3','\lambda_4','\lambda_5','location','sw')
% Axes handle 1 (this is the visible axes)
% ah1 = gca;
% Legend at axes 1
% lg1=legend(ah1,h(1:2),'$F^{\prime}_{Local}$ Param. 1','$F^{\prime}_{1}$ Param. 1','location','se');
% set(lg1,'interpreter','latex')
%
% Axes handle 2 (unvisible, only for place the second legend)
% ah2=axes('position',get(gca,'position'), 'visible','off');
% Legend at axes 2
% lg2=legend(ah2,h(3:7),'$\bar{\lambda}_1\,\,(x=0.1)$','$\bar{\lambda}_2\,\,(x=1.1)$','$\bar{\lambda}_3\,\,(x=2.1)$','$\bar{\lambda}_4\,\,(x=3.1)$','$\bar{\lambda}_5\,\,(x=4.1)$','location','sw');
lg2=legend('$\bar{\lambda}_1\,\,(x=0.1)$','$\bar{\lambda}_2\,\,(x=1.1)$','$\bar{\lambda}_3\,\,(x=2.1)$','$\bar{\lambda}_4\,\,(x=3.1)$','$\bar{\lambda}_5\,\,(x=4.1)$','location','sw');
set(lg2,'interpreter','latex')
%
figure()
plot(x(uu),Fprime1(uu),'bo');%,'MarkerSize',10);
hold on;
plot(x(uu),Fprime1_uu,'b+','MarkerSize',10);
plot(x(uu),Fprime2(uu),'rsq','MarkerSize',10);
plot(x(uu),Fprime2_uu,'rx','MarkerSize',10);
xlabel('x [m]');%,'interpreter','latex');
ylabel('Sesitivity Loads');
grid on;
lg1=legend('$F^{\prime}_{2}$ Param. 1','$F^{\prime}_{local}$ Param. 1','$F^{\prime}_{2}$ Param. 2','$F^{\prime}_{local}$ Param. 2','location','ne');
set(lg1,'interpreter','latex')
title('Sesitivity Loads');
%
figure()
plot(x(uu),Fprime1(uu),'bo');%,'MarkerSize',10);
hold on;
plot(x(uu),Fprime1_uu,'b+','MarkerSize',10);
plot(x(uu(2:end)),Fprime2(uu(2:end)),'rsq','MarkerSize',10);
plot(x(uu(2:end)),Fprime2_uu(2:end),'rx','MarkerSize',10);
xlabel('x [m]');%,'interpreter','latex');
ylabel('Sesitivity Loads');
grid on;
lg1=legend('$F^{\prime}_{2}$ Param. 1','$F^{\prime}_{local}$ Param. 1','$F^{\prime}_{2}$ Param. 2 (w/o first term)','$F^{\prime}_{local}$ Param. 2 (w/o first term)','location','ne');
set(lg1,'interpreter','latex')
title('Sesitivity Loads');