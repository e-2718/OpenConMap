clear all
clc

%%
alp=80;
beta=70;
F={@(x) tand(-alp).*x,@(x) tand(-alp)+0.*x,@(x) tand(alp-180).*x,@(x) tand(180-beta).*x,@(x) tand(180-beta)*-0.5+0.*x,@(x) tand(beta).*x};
X=[0,1+tand(-alp)*i,-1+tand(-alp)*i,0,-0.5+tand(beta)*0.5i,0.5+tand(beta)*0.5i];
z_V=Initialize(F,X);

figure
z=Boundary(linspace(0,1,5000));
plot(z,'color',[0, 114, 178]/255,'LineWidth',2);
hold on
plot(real([z_V,z_V(1)]),imag([z_V,z_V(1)]),'s-','color',[0, 158, 115]/255,'LineWidth',2,'MarkerSize',8);
axis image
axis([-5 5 -7 3])
grid on

h=gca;
h.FontSize=15;

%%
[C,err]=ConformalSolve(z_V,100,100,600,0.3,1e-4);

omegau=C(end:-1:1).';
omegad=[1 0].';

%%
[rho,theta]=meshgrid(linspace(0,1,10),linspace(0,2*pi,10000));
zeta1=rho.*exp(i*theta);
z1=(polyval(omegau,zeta1)./(polyval(omegad,zeta1)));

[rho,theta]=meshgrid(linspace(0,1,1000),linspace(0,2*pi,50));
zeta2=rho.*exp(i*theta);
z2=(polyval(omegau,zeta2)./(polyval(omegad,zeta2)));

zeta3=exp(i*linspace(0,2*pi,100));
z3=Boundary(linspace(0,1,5000));

H=figure;
H.Position=[10 10 200*2 155*2];
plot(z1,'color',[86, 180, 233]/255,'LineWidth',1.0);
hold on
plot(z2.','color',[86, 180, 233]/255,'LineWidth',1.0);
plot(z3,'color',[0, 114, 178]/255,'LineWidth',1.0);
axis image
axis([-5 5 -7 3])

h=gca;
h.FontSize=10;

H=figure;
H.Position=[10 10 200*2 150*2];
plot(zeta1,'color',[86, 180, 233]/255,'LineWidth',1.0);
hold on
plot(zeta2.','color',[86, 180, 233]/255,'LineWidth',1.0);
plot(zeta3,'color',[0, 114, 178]/255,'LineWidth',1.5);
axis image
axis([-1 1 -1 1])
h=gca;
h.FontSize=10;


%%
n1=20;
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

n1=20;
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

H=figure;
H.Position=[10 10 200*2 150*2];
plot(z1,'color',[230, 159, 0]/255,'LineWidth',1.0);
hold on
plot(z2,'color',[230, 159, 0]/255,'LineWidth',1.0);
plot(z3,'color',[213, 94, 0]/255,'LineWidth',1.0);
axis image
axis([-5 5 -7 3])

h=gca;
h.FontSize=10;

H=figure;
H.Position=[10 10 200*2 150*2];
plot(zeta1,'color',[230, 159, 0]/255,'LineWidth',1.0);
hold on
plot(zeta2,'color',[230, 159, 0]/255,'LineWidth',1.0);
plot(zeta3,'color',[213, 94, 0]/255,'LineWidth',1.5);
axis image
axis([-1 1 -1 1])
h=gca;
h.FontSize=10;














