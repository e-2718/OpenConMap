function [CC,epsilon]=ConformalSolve(zC,N,N2,N3,Lambda,er0)

Lambda1=0;
[Xi,etaB,C]=InitialValue_fun(zC,N2,N3,Lambda1);

gen=1;
epsilon=zeros(1,N+1);
epsilon(gen)=RelativeError(C);
epsilon0=1;
eps=epsilon(gen);

CC=C;
while and(epsilon0>=er0,gen<N+1)
    
    gen=gen+1
       
    zB=Boundary(etaB,'C');
    C=linsolve(Xi,zB);

    zA=Xi*C;

    rs=[0;cumsum(abs(zA(2:end)-zA(1:end-1)))];
    etaA=rs/(rs(end)+abs(zA(end)-zA(1)));

    etaB=etaA*Lambda+etaB*(1-Lambda);
    
%     if k1==N
    if 1
        epsilon(gen)=RelativeError(C);
        epsilon0=abs((epsilon(gen)-epsilon(gen-1))./epsilon(gen-1));
    end

    if epsilon(gen)<eps
        CC=C;
        eps=epsilon(gen);
    end
    
end

epsilon=epsilon(1:gen);

end






function [Xi,etaA,C,epslion0,theta]=InitialValue_fun(zC,N2,N4,Lambda1)

zC=[zC(2:end),zC(1)];
p=polygon(flip(zC));
f=extermap(p);

if isempty(f)
    
    [Xi,etaA]=fun1(N2,N4);
    C=[];
    epslion0=[];

elseif strcmp(class(f),'extermap')
    
    [Xi,etaA,C,epslion0,theta]=fun2(f,N2,N4,Lambda1);
    
end


end

function [Xi,etaA]=fun1(N2,N4)

    beta2=linspace(0,2*pi,N4+1);
    beta2(end)=[];
    
    Xi=exp(i*kron(beta2.',-1:N2));
    
    etaA=linspace(0,1,N4+1)';
    etaA(end)=[];
    
end

function [Xi,etaA,C,epslion0,theta]=fun2(f,N2,N4,Lambda1)

%%
% clf
% plot(r,t,'b.')
% hold on
% plot(etaB,beta1,'r.')
% 
% clf
% plot(abs(zA),'b.')
%%
    N3=19999;
    t=linspace(0,2*pi,N3);
    z=f(exp(1i*t));


    r0=abs(z(2:end)-z(1:end-1));
    r=[0,cumsum(r0)];
    r=r/r(end);

%     save Boundary_data3.mat r z

    %%
%     name3='Boundary_data3.mat';

    etaB=linspace(0,1,N4+1);
    etaB(end)=[];
%     zB=Boundary(etaB,name3);
%     XiB=evalinv(f,zB);

%     evalinv(f,zB(615));
%     evalinv(f,11.811141644046312 - 0.685203158017207*i);
%     
%     clf
%     plot(zB,'b.-')
%     hold on
%     plot(zB(615),'ro-')
    
    %%
%     beta1=angle(XiB);
    beta1=interp1(r*1e5,t,etaB*1e5);

%     plot(log10(r(2:end)-r(1:end-1)),'.')
% 
% sum(isnan(log10(etaB(2)-etaB(1))))
% 
% r(2)-r(1)
% 
% r(3334)
% r(3335)
% 
% t(10000)
% t(10001)
% 
% f(exp(1i*pi))
% f(exp(1i*t(10001)))
% %%
% clf
% plot(r,'.')
% hold on
% plot(r(1:end-1)-unique(r),'.')
    %%
    beta1(beta1<0)=beta1(beta1<0)+2*pi;
    beta2=linspace(0,2*pi,N4+1);
    beta2(end)=[];
    
    theta=beta1*Lambda1+(1-Lambda1)*beta2;

    Xi=exp(i*kron(theta.',-1:N2));
    C=Laurent_fun(f,N2);
    zA=Xi*C;

    etaA=[0;cumsum(abs(zA(2:end)-zA(1:end-1)))];
    etaA=etaA/(etaA(end)+abs(zA(end)-zA(1)));
    
    f_z=@(theta,Xi,C) Xi*C;
    epslion0=RelativeError(C);
%     epslion0=RelativeError(C);

    
end

function C = Laurent_fun(M,m)

p = polygon(M); 

beta = 1-angle(p);
r = parameters(M);
z = r.prevertex;
c = r.constant;

n = length(p);
gam = ones(n,m+2);
for k = 1:m+1
    gam(:,k+1) = -gam(:,k).*(beta-k+1)./(k*z);
end

e1 = zeros(m+2,1);
e1(1) = 1;
C = -c*e1;
for j = 1:n
    C = toeplitz(gam(j,:).', e1')*C;
end
C = C(3:m+2)./(-(1:m)');
x0 = 10^(-10/(m+1));
c0 = eval(M,x0) + c/x0 - x0.^(1:m)*C;
C = [c0;C];
C=[-c;C];

end








