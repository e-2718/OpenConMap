function za=BoundaryInitial(F,X0,N5,H,P)
N1=length(X0);
H=0;
P='A';

ss=[];
for k1=1:N1
    fs=func2str(F{k1});
    ss=[ss,fs(3)];
end

X=[X0,X0(1)];

%%
z=[];
zz=cell(1,N1);
d=0.001;

for k1=1:N1

    A=X(k1);
    B=X(k1+1);

    N3=round(abs(A-B)/d);
    N2=N3+1;
    while N3~=N2
        N2=N3;
            if ss(k1)=='x'
                
                x=linspace(real(A),real(B),N2);
                y=x+i*F{k1}(x);
                N3=round(sum(abs(y(2:end)-y(1:end-1)))/d);
                
            elseif ss(k1)=='y'
                
                x=linspace(imag(A),imag(B),N2);
                y=F{k1}(x)+i*x;
                N3=round(sum(abs(y(2:end)-y(1:end-1)))/d);
                
            end
    end

    a=abs(y(2:end-1)-y(1:end-2));
    b=abs(y(3:end)-y(2:end-1));
    c=abs(y(3:end)-y(1:end-2));
    p=(a+b+c)/2;
    s=real((p.*(p-a).*(p-b).*(p-c)).^0.5);
    R{k1}=(4.*s)./(a.*b.*c);
    r(k1)=sum(R{k1});
    zz{k1}=y;
    
    z=[z,y(1:end-1)];

end

z=[z,z(1)];

%%
NN1=round(r./sum(r)*(N5-N1));
NN2=r./(NN1+1);
za=[];

for k1=1:N1
    
    za=[za,X(k1)];
    
    r0=linspace(0,1,NN1(k1)+2)*r(k1);
    rr=[0,cumsum(R{k1})];

    for k2=2:length(r0)-1
        [~,u]=min(abs(rr-r0(k2)));
        za=[za,zz{k1}(u)];
    end

end
%%
al=angle(conj(z));

if P=='B'
    al(al<pi/2)=al(al<pi/2)+2*pi;
    al(1)=5*pi/2;
elseif P=='A'
    al(al<0)=al(al<0)+2*pi;
    al(end)=2*pi;
end

% plot(al*180/pi,'.-k')
% 
% clf 
% plot(real(z),imag(z),'.-k')
% axis image
% hold on
% axis([-6 6 [-4 4]-H])
% for k1=round(linspace(1,length(z),100))
%     plot(z(k1),'*k')
%     pause(0.003)
% end

save Boundary_data2.mat al z H
%%
r0=abs(z(2:end)-z(1:end-1));
r=[0,cumsum(r0)];
r=r/r(end);

save Boundary_data1.mat r z H
%%
end































