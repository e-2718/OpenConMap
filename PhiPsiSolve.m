function [phiu,phid,psiu,psid,omegau,omegad,a] = PhiPsiSolve(C,g1,g2,P,Q,Alpha)
C=C.';
g2=g2(end:-1:1);
N=size(C,2)-2;

%%
n1=C(end:-1:1);
n2=1:N;

W=zeros(N,2*N+1);
for k1=1:N
    n3=N-k1;
    n4=zeros(1,n3+1);
    n4(1)=k1;
    n6=conv(n1,n4);
    W(k1,1:size(n6,2))=n6(end:-1:1);
end
W=W(:,end:-1:1);

U=[N+2,17,N-[1:N]+1];
U=[-conj(C(1)),0,conj(C(3:N+2)).*[1:N],0];

D=[];
for k1=1:N
    D=[D;deconv(W(k1,:),U)];
end
D=D(:,end:-1:1);

D=[D(1,:)*0;D];
D=[D,zeros(N+1,2)];
%%
D1=real(D);
D2=imag(D);
D3=-real(D);

for k1=1:N+1
    D1(k1,k1)=D1(k1,k1)+1;
end
for k1=1:N+1
    D3(k1,k1)=D3(k1,k1)+1;
end

G=[D1,D2;D2,D3];

BB=inv(G.')*[real(g1);imag(g1)];

a=BB(1:N+1)+BB(N+2:N*2+2)*i;
phi=(Q+P)/4*C.';
phi(2:end)=phi(2:end)+a;

phiu=phi(end:-1:1);
phid=[1 0].';
%%
n1=conj(C);
n2=[1:N].*a(2:end).';
n2=n2(end:-1:1);

W=conv(n1,n2);
U=[-C(1),0,C(3:N+2).*[1:N]];
U=U(end:-1:1);
U=[U,zeros(1,N-2)];

[m1,m2] = deconv(W,U);
U1=U(end:-1:1);
M=m2(end:-1:1);

U1(1:N-2)=[];
M(end-1:end)=[];

[m3,m4] = deconv(M,U1);
m4=m4(end:-1:1);
m4=m4(1:N+1);
U1=U1(end:-1:1);

%%
psi0=g2.';
psi0(end)=psi0(end)-conj(a(1));
psi0(end-1:end)=psi0(end-1:end)-m1;
psi1=-1/2*(P - Q)*exp(-2*i*Alpha)*C;
psi1=psi1(end:-1:1);

psi1(1:N+1)=psi1(1:N+1).'+psi0.';
psi1=psi1(end-1:end);
psiu=-m4;
psid=U1;

s1=conv(psi1,psid);
s1(2:end)=s1(2:end)+[psiu 0];
psiu=s1.'/psid(1);
psid=[psid 0].'/psid(1);
%%
omegau=C(end:-1:1).';
omegad=[1 0].';
a=a(end:-1:1);
end

