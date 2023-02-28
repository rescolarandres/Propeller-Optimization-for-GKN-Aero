function [Ct,Cp,eta,T,P,V] = BET_RE(j,chord,r,beta,mission,blade)
%This function computes thrust and power using the BET
%The inputs are blade geometry and working conditions
% For fixed omega
    omega = mission.omega;
    n = mission.omega/(2*pi);
    V = blade.diam*j*n;
%For fixed V
%     V = mission.V;
%     n = V/(blade.diam*j);
%     omega = n*2*pi;
    gamma=atan(V./(r*omega));
    alpha = beta-gamma;
    for k=1:length(r)
        Re(k) = mission.rho*chord(k)*V/mission.viscosity;
        [Cl(k),Cd(k)]=get_coeffs(rad2deg(alpha(k)),Re(k),blade);
%         [Cl,Cd]=Aer_Coeff_NACA0012('linear',alpha(k)*180/pi,Re(k));
    end
    dLdr = 0.5*mission.rho.*Cl.*chord.*(V^2+omega^2*r.^2);
    dDdr = 0.5*mission.rho.*Cd.*chord.*(V^2+omega^2*r.^2);
    dTdr = blade.num*(dLdr.*cos(gamma)-dDdr.*sin(gamma));
    dPdr = omega*r*blade.num.*(dLdr.*sin(gamma)+dDdr.*cos(gamma));
    T = trapz(r,dTdr);
    P = trapz(r,dPdr);
    Ct = T/(mission.rho*blade.diam^4*n^2);
    Cp = P/(mission.rho*blade.diam^5*n^3);
    eta = T*V./P;
    
    
end