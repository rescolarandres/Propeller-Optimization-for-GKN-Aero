function [Ct,Cp,eta,T,P,Ve,dTdr,dPdr,alpha,cl,cd]=BEMT_RE(j,mission,blade,r,beta,omega,rho)
% This function calculates Thrust T[N] and Power P [W] and their
% derivatives along the blade span (respectively dTdr [N/m] and dPdr [W/m])
% using the Blade-Element Momentum Theory.

% Given rpms
n = omega/(2*pi);
V = blade.diam*j*n;

r = r*blade.tip_rad;
% Finding the gamma distribution for the given working conditions
options = optimset('TolFun',1e-6);
gamma=fsolve(@(gamma) induction(gamma,V,blade,beta,blade.chord,r,mission,omega,rho),beta,options);

%Calculating Thrust and Power distributions given the correct gamma
[~,dTdr,dPdr,Ve,alpha,cl,cd,a,b]=induction(gamma,V,blade,beta,blade.chord,r,mission,omega,rho);

%Integration of Thrust and Power along the blade
T=trapz(r,dTdr);
P=trapz(r,dPdr);
Ct = T/(rho*blade.diam^4*n^2);
Cp = P/(rho*blade.diam^5*n^3);
eta = T*V./P;

end

function [error,dTdr,dPdr,Ve,alpha,cl,cd,a,b]=induction(gamma,V,blade,beta,chord,r,mission,omega,rho)
% This function solves for the correct gamma distribution at given working
% conditions, and calculates Thrust and Power distributions (dTdr and
% dPdr).

sol=blade.num.*chord./(2.*pi.*r); % propeller solidity distribution
alpha=beta-gamma; % angle of attack distribution

Ve=(V.^2+(omega*r).^2).^0.5; % resultant element velocity distribution

[~,lambda1,lambda2,a,b,cl,cd]=reynolds(Ve,blade,gamma,sol,r,chord,mission,alpha,omega,rho,V);

% Thrust and Power distributions
dTdr=0.5.*rho.*blade.num.*Ve.^2.*lambda1.*chord;
dPdr=0.5.*rho.*blade.num*omega.*Ve.^2.*lambda2.*chord.*r;

%Error function 1 (correct gamma is a root of it)
error=gamma-atan(V./(omega.*r).*(1+a)./(1-b));
end

function [Ve,lambda1,lambda2,a,b,Cl,Cd]=reynolds(Ve,blade,gamma,sol,r,chord,mission,alpha,omega,rho,V)
% This function solves for the correct Reynolds number (including propeller
% auto-induction effect). 
% It also calculates:
% - Ve: effective velocity
% - lambda1: projection of aerodynamic forces coefficients in the axis of
%   the propeller.
% - lambda2: projection of aerodynamic forces coefficients in the plane of
%   the propeller.
% - a: axial induction ratio (a=v/V)
% - b: rotational induction ratio (b=w/omega)
% Working fluid is considered to be air.
Re = rho*Ve.*chord/mission.viscosity;
Cl=0*r;Cd=Cl;
for k=1:length(r)
        [Cl(k),Cd(k)]=get_coeffs(rad2deg(alpha(k)),Re(k),blade);
end

lambda1=Cl.*cos(gamma)-Cd.*sin(gamma);
lambda2=Cl.*sin(gamma)+Cd.*cos(gamma);

subs = sol./(4*sin(gamma).^2).*lambda1;
a=subs./(1-subs);

subs = sol./(4*sin(gamma).*cos(gamma)).*lambda2;
b=subs./(1+subs);

Ve=(V.^2.*(1+a).^2+omega.^2.*r.^2.*(1-b).^2).^0.5;


end

