clear all
clc

%% Draw physical domain boundarie
y1=linspace(1,-1,100);
f=@(y) 0.5*(1-y.^2).^0.5;
x1=f(y1);

x2=linspace(0,-1,100);
f=@(x) -1+x.*0;
y2=f(x2);

y3=linspace(-1,0,100);
f=@(y) -1+y.*0;
x3=f(y3);

x4=linspace(-1,0,100);
f=@(x) x+1;
y4=f(x4);

figure
hold on
plot(x1,y1,'-','color',[0, 158, 115]/255,'LineWidth',2,'MarkerSize',8);
plot(x2,y2,'-','color',[0, 158, 115]/255,'LineWidth',2,'MarkerSize',8);
plot(x3,y3,'-','color',[0, 158, 115]/255,'LineWidth',2,'MarkerSize',8);
plot(x4,y4,'-','color',[0, 158, 115]/255,'LineWidth',2,'MarkerSize',8);

axis image

%% Image of polygon approximation 
F={@(y) 0.5*(1-y.^2).^0.5, @(x) -1, @(y) -1, @(x) x+1};
X=[i, -i, -1-i, -1];

z_V=Initialize(F, X, 10);

figure
z=Boundary(linspace(0,1,5000));
plot(z,'color',[0, 114, 178]/255,'LineWidth',2);
hold on
plot(real([z_V,z_V(1)]),imag([z_V,z_V(1)]),'s-','color',[0, 158, 115]/255,'LineWidth',2,'MarkerSize',8);
axis image
axis([-1.2 0.5+0.2 -1.2 1.2])
grid on

%% Solve the mapping function
[C, err] = ConformalSolve(z_V, 50, 30, 200, 0.6, 1e-3);

omegau=C(end:-1:1).';
omegad=[1 0].';

%% Image of the mapping function
[rho,theta]=meshgrid(linspace(1,5,50),linspace(0,2*pi,300));
rho=1./rho;
zeta1=rho.*exp(i*theta);
z=(polyval(omegau,zeta1)./(polyval(omegad,zeta1)));

mesh(real(z),imag(z),real(z).*0)
view([0 0 1])
axis image

%% Solve the stress and displacement
E=2e9; % 弹性模量
niu=0.2; % 泊松比
G=E/(2*niu+2); % 剪切模量

sx=15; % 水平应力
sy=20; % 垂直应力
sxy=3; % 剪切应力

[P,Q,Alpha]=StressCondition(sx,sy,sxy);
[g1,g2]=EquivalentCondition(P,Q,Alpha,C);
[phiu,phid,psiu,psid,omegau,omegad] = PhiPsiSolve(C,g1,g2,P,Q,Alpha);
[sigma_x,sigma_y,tau_xy,u_x,u_y] = SDSolve(phiu,phid,psiu,psid,omegau,omegad,rho,theta,G,niu);

%% Image of the stress
mesh(real(z),imag(z),sigma_x)















