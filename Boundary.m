function zz=Boundary(r0,R,name)

if nargin==1
    R='C';
end

if R=='C';
    
    if nargin<3
        name='Boundary_data1.mat';
    end
    load(name)
    zz=interp1(r',z.',r0,'linear')-H*i;

elseif R=='P';
    
    if nargin<3
        name='Boundary_data2.mat';
    end
    load(name)
    zz=interp1(al.',z.',r0,'linear')-H*i;
    
end
% 
% 
% if nargin==1
%     name='Boundary_data1.mat';
% end
% 
% load(name)
% zz=interp1(r',z.',r0,'linear')-H*i;

end

% plot(al.','b.-')
% plot(r0,'r.-')

% plot(abs(zz(2:end)-z(1:end-1)),'b.-')









