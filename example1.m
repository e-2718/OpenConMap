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

% F={@(y) 0.5*(1-y.^2).^0.5, @(x) -1, @(y) -1, @(y) y-1};
% X=[i, -i, -1-i, -1];

% F={@(x) -1, @(y) -1, @(x) x+1, @(y) 0.5*(1-y.^2).^0.5};
% X=[-i, -1-i, -1, i];

% F={@(y) -1, @(x) x+1, @(y) 0.5*(1-y.^2).^0.5, @(x) -1};
% X=[-1-i, -1, i, -i];

% F={@(x) x+1, @(y) 0.5*(1-y.^2).^0.5, @(x) -1, @(y) -1};
% X=[-1, i, -i, -1-i];

z_V=Initialize(F, X, 10);

figure
z=Boundary(linspace(0,1,5000));
plot(z,'color',[0, 114, 178]/255,'LineWidth',2);
hold on
plot(real([z_V,z_V(1)]),imag([z_V,z_V(1)]),'s-','color',[0, 158, 115]/255,'LineWidth',2,'MarkerSize',8);
axis image
axis([-1.2 0.5+0.2 -1.2 1.2])
grid on

%% solve
[C, err] = ConformalSolve(z_V, 50, 30, 200, 0.6, 1e-3);

omegau=C(end:-1:1).';
omegad=[1 0].';

%%
figure
plot(err,'.-','color',[0, 114, 178]/255,'LineWidth',1.5);
hold on
grid on
h=gca;
h.YScale='log';
h.FontSize=15;
h.XAxis.Color=[0, 114, 178]/255;
h.YAxis.Color=[0, 114, 178]/255;

%% The image of an orthogonal grid under conformal mapping.
[rho,theta]=meshgrid(linspace(0,1,20),linspace(0,2*pi,10000));
zeta1=rho.*exp(i*theta);
z1=(polyval(omegau,zeta1)./(polyval(omegad,zeta1)));

[rho,theta]=meshgrid(linspace(0,1,100),linspace(0,2*pi,100));
zeta2=rho.*exp(i*theta);
z2=(polyval(omegau,zeta2)./(polyval(omegad,zeta2)));

zeta3=exp(i*linspace(0,2*pi,100));
z3=Boundary(linspace(0,1,5000));

figure
plot(zeta1,'color',[86, 180, 233]/255,'LineWidth',1.0);
hold on
plot(zeta2.','color',[86, 180, 233]/255,'LineWidth',1.0);
plot(zeta3,'color',[0, 114, 178]/255,'LineWidth',2.0);
axis image

figure
plot(z1,'color',[86, 180, 233]/255,'LineWidth',1.0);
hold on
plot(z2.','color',[86, 180, 233]/255,'LineWidth',1.0);
plot(z3,'color',[0, 114, 178]/255,'LineWidth',1.5);

axis image
axis([-2 2 -2 2])

%% The image of another orthogonal grid under conformal mapping.
n1=40;
n2=2000;
xt0=linspace(-1,1,n1);
xt1=[];
yt1=[];
for k1=1:n1
    yt1=[yt1;linspace(-sqrt(1-xt0(k1).^2),sqrt(1-xt0(k1).^2),n2)];
end
yt1=yt1';
for k1=1:n2
    xt1=[xt1;xt0];
end
[theta,rho] = cart2pol(xt1,yt1);
zeta1=rho.*exp(i*theta);
z1=(polyval(omegau,zeta1)./(polyval(omegad,zeta1)));

n1=40;
n2=2000;
yt0=linspace(-1,1,n1);
xt2=[];
yt2=[];
for k1=1:n1
    xt2=[xt2;linspace(-sqrt(1-yt0(k1).^2),sqrt(1-yt0(k1).^2),n2)];
end
xt2=xt2';
for k1=1:n2
    yt2=[yt2;yt0];
end
[theta,rho] = cart2pol(xt2,yt2);
zeta2=rho.*exp(i*theta);
z2=(polyval(omegau,zeta2)./(polyval(omegad,zeta2)));

zeta3=exp(i*linspace(0,2*pi,100));
z3=Boundary(linspace(0,1,5000));

figure
plot(zeta1,'color',[230, 159, 0]/255,'LineWidth',1.0);
hold on
plot(zeta2,'color',[230, 159, 0]/255,'LineWidth',1.0);
plot(zeta3,'color',[213, 94, 0]/255,'LineWidth',2.0);
axis image

figure
plot(z1,'color',[230, 159, 0]/255,'LineWidth',1.0);
hold on
plot(z2,'color',[230, 159, 0]/255,'LineWidth',1.0);
plot(z3,'color',[213, 94, 0]/255,'LineWidth',1.5);

axis image
axis([-2 2 -2 2])




