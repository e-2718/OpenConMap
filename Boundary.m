function z=Boundary(x,R,name)

if nargin==1
    R='C';
end

if R=='C';
    
    if nargin<3
        name='Boundary_data1.mat';
    end
    load(name)
    z=interp1(r',z.',x,'linear')-H*i;

elseif R=='P';
    
    if nargin<3
        name='Boundary_data2.mat';
    end
    load(name)
    z=interp1(al.',z.',x,'linear')-H*i;
    
end

end









