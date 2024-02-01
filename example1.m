clear all
clc

%%
F={@(y) 0.5*(1-y.^2).^0.5, @(x) -1, @(y) -1, @(x) x+1};
X=[i, -i, -1-i, -1];

% F={@(y) 0.5*(1-y.^2).^0.5, @(x) -1, @(y) -1, @(y) y-1};
% X=[i, -i, -1-i, -1];

% F={@(x) -1, @(y) -1, @(x) x+1, @(y) 0.5*(1-y.^2).^0.5};
% X=[-i, -1-i, -1, i];

% F={@(y) -1, @(x) x+1, @(y) 0.5*(1-y.^2).^0.5, @(x) -1};
% X=[-1-i, -1, i, -i];

% F={@(x) x+1, @(y) 0.5*(1-y.^2).^0.5, @(x) -1, @(y) -1};
% X=[-1, i, -i, -1-i];

z = Initialize(F, X, 10);
[C, err] = ConformalSolve(z, 50, 30, 200, 0.6, 1e-3);

%%
clf
zz=Boundary(linspace(0,1,5000));
plot(real(zz),imag(zz),'color',[0, 114, 178]/255,'LineWidth',2);
hold on
plot(real([z,z(1)]),imag([z,z(1)]),'s-','color',[0, 158, 115]/255,'LineWidth',2,'MarkerSize',8);
axis image
axis([-1.2 0.5+0.2 -1.2 1.2])
grid on

h=gca;
h.FontSize=15;
h.FontName='微软雅黑';
h.XAxis.Color=[0, 114, 178]/255;
h.YAxis.Color=[0, 114, 178]/255;
%%
omegau=C(end:-1:1).';
omegad=[1 0].';
%
[rho,theta]=meshgrid(linspace(1,3,20),linspace(0,2*pi,10000));
rho=1./rho;
zeta=rho.*exp(i*theta);
z=(polyval(omegau,zeta)./(polyval(omegad,zeta)));

x=real(z);
y=imag(z);

clf
plot(x,y,'color',[213, 94, 0]/255,'LineWidth',1.5);
hold on

[rho,theta]=meshgrid(linspace(1,3,100),linspace(0,2*pi,100));
rho=1./rho;
zeta=rho.*exp(i*theta);
z=(polyval(omegau,zeta)./(polyval(omegad,zeta)));

x=real(z);
y=imag(z);

plot(x',y','color',[213, 94, 0]/255,'LineWidth',1.5);
axis image

zz=Boundary(linspace(0,1,5000));
plot(real(zz),imag(zz),'color',[0, 114, 178]/255,'LineWidth',1.5);

N=0.8;
axis([-1.-N 0.5+N -1.-N 1.0+N])

h=gca
h.FontSize=15
h.FontName='微软雅黑'
h.XAxis.Color=[0, 114, 178]/255;
h.YAxis.Color=[0, 114, 178]/255;
%%
clf
plot(err,'s-','color',[0, 114, 178]/255,'LineWidth',1.5);
hold on
grid on
h=gca;
h.YScale='log';
h.FontSize=15;
h.FontName='微软雅黑';
h.XAxis.Color=[0, 114, 178]/255;
h.YAxis.Color=[0, 114, 178]/255;




