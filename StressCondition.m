function [p0,q0,alpha0]=StressCondition(sx0,sy0,sxy0)
%%
    if sxy0==0
        p0=sx0;
        q0=sy0;
        alpha0=0;
    elseif sx0==sy0
        p0=(sx0+sy0)/2+(((sx0-sy0)/2)^2+sxy0^2)^0.5;
        q0=(sx0+sy0)/2-(((sx0-sy0)/2)^2+sxy0^2)^0.5;
        alpha0=atan((-sx0+sy0)/(2*sxy0)+(((sx0-sy0)/(2*sxy0))^2+1)^0.5);
        if sxy0<0
            alpha0=alpha0-pi/2;
        end
    else
        p0=(sx0+sy0)/2+sign(sx0-sy0)*(((sx0-sy0)/2)^2+sxy0^2)^0.5;
        q0=(sx0+sy0)/2-sign(sx0-sy0)*(((sx0-sy0)/2)^2+sxy0^2)^0.5;
        alpha0=atan((-sx0+sy0)/(2*sxy0)+(((sx0-sy0)/(2*sxy0))^2+1)^0.5);
        if sx0<sy0
            if sxy0>0
                alpha0=alpha0-pi/2;
            end
        else
            if sxy0<0
                alpha0=alpha0-pi/2;
            end
        end
    end

end