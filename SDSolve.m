function [sigma_x,sigma_y,tau_xy,u_x,u_y] = SDSolve(phiu,phid,psiu,psid,omegau,omegad,rho,theta,G,niu)
C=omegau(end:-1:1).';
N=size(C,2)-2;
zeta=rho.*exp(i*theta);
% z=(polyval(omegau,zeta)./(polyval(omegad,zeta)));

[q1,q2]=polyder(phiu,phid);
[q3,q4]=polyder(omegau,omegad);
[q5,q6]=polyder(psiu,psid);
[q7,q8]=polyder(q1,q3);

omega_data=polyval(omegau,zeta)./(polyval(omegad,zeta));
omega_conj=conj(omega_data);
omega_diff=polyval(q3,zeta)./(polyval(q4,zeta));
omega_conj_diff=conj(omega_diff);

phi_data=polyval(phiu,zeta)./(polyval(phid,zeta));
phi_diff=polyval(q1,zeta)./(polyval(q2,zeta));
phi_conj_diff=conj(phi_diff);

psi_diff=polyval(q5,zeta)./(polyval(q6,zeta));
psi_conj=conj(polyval(psiu,zeta)./(polyval(psid,zeta)));

phi_diff_psi_diff_diff=polyval(q7,zeta)./(polyval(q8,zeta));

p1=4*real(phi_diff./omega_diff);
p2=2*zeta.^2./(rho.^2.*omega_conj_diff).*(omega_conj.*phi_diff_psi_diff_diff+psi_diff);
p3=-exp(-i*theta).*omega_conj_diff./(2*G*abs(omega_diff)).*((3-4*niu).*phi_data-omega_data./omega_conj_diff.*phi_conj_diff-psi_conj);

% p1=4*real((polyval(q1,zeta))./(polyval(q3,zeta)));
% p2=2*zeta.^2./(abs(zeta).^2.*conj(polyval(q3,zeta)./(polyval(q4,zeta)))).*(conj(polyval(omegau,zeta)./(polyval(omegad,zeta))).*(polyval(q7,zeta)./(polyval(q8,zeta)))+(polyval(q5,zeta)./(polyval(q6,zeta))));

tau=imag(p2)/2;
sigma_rho=(p1-real(p2))/2;
sigma_theta=(p1+real(p2))/2;

u_rho=real(p3);
u_theta=imag(p3);

h1=-rho.^(-2).*(real(C(1)).*cos(theta)+imag(C(1))*sin(theta));
for k1=1:N
    h1=h1+k1*rho.^(k1-1).*(real(C(k1+2)).*cos(k1*theta)-imag(C(k1+2)).*sin(k1*theta));
end

h2=0;
for k1=-1:N
    h2=h2+rho.^k1.*(-k1*real(C(k1+2)).*sin(k1*theta)-k1*imag(C(k1+2)).*cos(k1*theta));
end

h3=-rho.^(-2).*(-real(C(1)).*sin(theta)+imag(C(1))*cos(theta));
for k1=1:N
    h3=h3+k1*rho.^(k1-1).*(real(C(k1+2)).*sin(k1*theta)+imag(C(k1+2)).*cos(k1*theta));
end

h4=0;
for k1=-1:N
    h4=h4+rho.^k1.*(k1*real(C(k1+2)).*cos(k1*theta)-k1*imag(C(k1+2)).*sin(k1*theta));
end

aa11=h1./(h1.^2+h3.^2).^0.5;
aa12=h3./(h1.^2+h3.^2).^0.5;

aa21=h2./(h2.^2+h4.^2).^0.5;
aa22=h4./(h2.^2+h4.^2).^0.5;

sigma_x=aa11.*(aa11.*sigma_rho + aa21.*tau) + aa21.*(aa11.*tau + aa21.*sigma_theta);
sigma_y=aa12.*(aa12.*sigma_rho + aa22.*tau) + aa22.*(aa12.*tau + aa22.*sigma_theta);
tau_xy=aa12.*(aa11.*sigma_rho + aa21.*tau) + aa22.*(aa11.*tau + aa21.*sigma_theta);

u_x=aa11.*u_rho+aa12.*u_theta;
u_y=aa21.*u_rho+aa22.*u_theta;

end

