function [g1,g2]=EquivalentCondition(P,Q,Alpha,C)

g1=[-1/2*((Q+P)*C(2)+(Q-P)*conj(C(2))*exp(2*i*Alpha)) 
        -1/2*((Q+P)*C(3)+(Q-P)*conj(C(1))*exp(2*i*Alpha))
        -1/2*(Q+P)*C(4:end)];

g2=[-1/2*((Q+P)*conj(C(2))+(Q-P)*C(2)*exp(-2*i*Alpha)) 
        -1/2*((Q+P)*conj(C(1))+(Q-P)*C(3)*exp(-2*i*Alpha))
        -1/2*(Q-P)*C(4:end)*exp(-2*i*Alpha)];
    
end

