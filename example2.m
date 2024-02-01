clear all
clc

%%
F={@(x) sin(x)+0.5, @(x) -x+6*pi, @(x) sin(x)-0.5, @(x) -x-6*pi};
X=[-19.100874546062855 + 0.251318624524098i
        18.598237297014663 + 0.251318624524097i
        19.100874546062855 - 0.251318624524098i
       -18.598237297014663 - 0.251318624524096i]';
z=Initialize(F,X,80);
[C,err]=ConformalSolve(z,1000,500,4000,0.6,1e-4);

X=[-19.100874546062855 + 0.251318624524098i 18.598237297014663 + 0.251318624524097i 19.100874546062855 - 0.251318624524098i -18.598237297014663 - 0.251318624524096i]';
%%
clf
zz=Boundary(linspace(0,1,5000));
plot(real(zz),imag(zz),'color',[0, 114, 178]/255,'LineWidth',2);
hold on
plot(real([z,z(1)]),imag([z,z(1)]),'s-','color',[0, 158, 115]/255,'LineWidth',2,'MarkerSize',8);
axis image
axis([-20 20 -2 2])
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

axis([-25 25 -5 5])

h=gca;
h.FontSize=15;
h.FontName='微软雅黑';
h.XAxis.Color=[0, 114, 178]/255;
h.YAxis.Color=[0, 114, 178]/255;
%%
clf
plot(err,'.-','color',[0, 114, 178]/255,'LineWidth',1.5);
hold on
grid on
h=gca;
h.YScale='log';
h.FontSize=15;
h.FontName='微软雅黑';
h.XAxis.Color=[0, 114, 178]/255;
h.YAxis.Color=[0, 114, 178]/255;




