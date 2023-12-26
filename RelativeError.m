function [e,ee,eta]=RelativeError(F,Z)
%%
theta=[linspace(0,2*pi,20000)].';
z=F(theta,Z{1},Z{2});

r=abs(z(2:end)-z(1:end-1))';
r=[0,cumsum(r)];
r=r/r(end);

H=0;
save Boundary_data5.mat r z H

%%
eta=linspace(0,1,5000);
z1=Boundary(eta,'C','Boundary_data5.mat');
z2=Boundary(eta,'C','Boundary_data1.mat');

%%
% zz=z2;
% clf 
% plot(zz,'.k')
% axis image
% hold on
% 
% axis([-1 1 -1 1]*1*1.5)
% for k1=round(linspace(1,length(zz),100))
%     plot(zz(k1),'*b')
%     pause(0.01)
% end


%%
x1=real(z1(1:end-1));
x2=real(z1(2:end));
x3=real(z2(2:end));
x4=real(z2(1:end-1));

y1=imag(z1(1:end-1));
y2=imag(z1(2:end));
y3=imag(z2(2:end));
y4=imag(z2(1:end-1));

l1=abs((x2-x1)+i*(y2-y1));
% l2=abs((x4-x3)+i*(y4-y3));

ee=abs((x1-x3).*(y2-y4)-(y1-y3).*(x2-x4))/2;
e=sum(ee)/sum(l1);
ee=ee./l1/sum(l1);

end




