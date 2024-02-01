clear all
clc

%%
C=[-5.400/2 -sind(87)*3.656/2];
D=[5.400/2 -sind(87)*3.656/2];

A=[cosd(180-87),sind(180-87)]*2.952+D;
B=[cosd(87),sind(87)]*3.656+C;

xx=[A(1),B(1),C(1),D(1)];
yy=[A(2),B(2),C(2),D(2)];

clf
plot(xx,yy,'s-')
hold on
axis image

pAB=polyfit([A(1),B(1)],[A(2),B(2)],1);
pBC=polyfit([C(1),B(1)],[C(2),B(2)],1);
pCD=polyfit([C(1),D(1)],[C(2),D(2)],1);
pDA=polyfit([A(1),D(1)],[A(2),D(2)],1);

%%
F0={@(x) pCD(1).*x+pCD(2), @(x) pBC(1).*x+pBC(2), @(x) pAB(1).*x+pAB(2),@(x) pDA(1).*x+pDA(2)};
xx=[D(1),C(1),B(1),A(1)];
yy=[D(2),C(2),B(2),A(2)];
X0=xx+yy*i;

clf
plot(xx,yy,'s')
hold on
axis image

xx=[xx,xx(1)];
for k1=1:length(X0)
    xs=linspace(xx(k1),xx(k1+1),100);
    ys=F0{k1}(xs);
    plot(xs,ys,'.')
end

%%
N1=50; % Iteration time
N2=50; % Series order
N3=max(N2*8,500); % Number of corresponding points
N4=4; % Number of vertices
Lambda=0.6; % Iteration speed
eff=1e-3; % Efficiency

zC=Initialize(F0,X0,N4);

eta=linspace(0,1,5000);
zz=Boundary(eta);

clf 
plot(zz,'k')
axis image
hold on

for k1=round(linspace(1,5000,100))
    plot(zz(k1),'*k')
    pause(0.01)
end

plot(zC,'sr-','LineWidth',2)
[C,err]=ConformalSolve(zC,N1,N2,N3,Lambda,eff);

%%
H1=figure;
H1.Position=[10 10 400 400];
H2=axes;
H2.Position=[0,0,1,1];

omegau=C(end:-1:1).';
omegad=[1 0].';

ss=2;
tt=2.5;
W=8;

% color1=[0, 114, 178]/255;
color1=[213, 94, 0]/255;
% color1=[86, 180, 233]/255;

[rho,theta]=meshgrid(linspace(1,ss,4),linspace(0,2*pi,500));
rho=1./rho;
zeta=rho.*exp(i*theta);
z=(polyval(omegau,zeta)./(polyval(omegad,zeta)));

x1=real(z);
y1=imag(z);

[rho,theta]=meshgrid(linspace(1,ss,200),linspace(0,2*pi,13));
rho=1./rho;
zeta=rho.*exp(i*theta);
z=(polyval(omegau,zeta)./(polyval(omegad,zeta)));

x2=real(z);
y2=imag(z);

[rho,theta]=meshgrid(linspace(ss,tt,200),linspace(0,2*pi,13));
rho=1./rho;
zeta=rho.*exp(i*theta);
z=(polyval(omegau,zeta)./(polyval(omegad,zeta)));

x3=real(z);
y3=imag(z);

axis([-1 1 -1 1]*6.7)

h=fill(x1(:,end),y1(:,end),color1);
h.EdgeColor='none';

hold on

h=fill(x1(:,1),y1(:,1),'w');
h.EdgeColor='none';

plot(x1(:,1:end-1),y1(:,1:end-1),'w','LineWidth',W);
plot(x2',y2','w','LineWidth',W);
plot(x3',y3','color',color1,'LineWidth',W);

view([0 0 1])
axis image

axis([-1 1 -1 1]*6.7)

axis off






