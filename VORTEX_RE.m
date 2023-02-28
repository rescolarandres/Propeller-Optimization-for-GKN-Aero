function [Ct,Cp,eta,T,P,Vef_nz,dTdr,dQdr,w] = VORTEX_RE(j,mission,blade,r,beta,omega,rho)
% This function calculates Thrust T[N] and Power P [W] and their
% derivatives along the blade span (respectively dTdr [N/m] and dPdr [W/m])
% using the Vortex Theory.

n = omega/(2*pi);
V = blade.diam*j*n;

% Finding the Circulation/AoA distribution for the given working conditions
% using Tremmel method
tol = 1e-4; % Tolerance set for the relaxation algorithm

alpha_i0= 0.04*r;   % Initial condition of induced AOA
circ0 = [0 10 14 18*ones(1,length(r)-6) 14 10 0];    % Initial condition of Circulation
error_circ = 1;
error_alpha=1;
n_iter = 0;
while error_alpha > tol && n_iter < 8e3
    while error_circ > tol
        % First, the circulation has to converge to then iterate the AoA 
        [~,~,~,~,~,~,circ] = circulation_sweep(alpha_i0,circ0,V,r,blade,mission,beta,omega,rho);
        circ0 = circ0 + 0.2*(circ-circ0);   % Under relaxation factor of 0.2 is applied
        error_circ = norm(circ-circ0);
    end
    % Secondly, the AoA is iterated until it converges, therefore obtaining the final solution
    [Vef_nz,alpha_ef,cl,cd,w,alpha_i,~] = circulation_sweep(alpha_i0,circ0,V,r,blade,mission,beta,omega,rho);
    alpha_i0 = alpha_i0 + 0.8*(alpha_i-alpha_i0);
    error_alpha = norm(alpha_i-alpha_i0);
    n_iter = n_iter + 1;
end
if n_iter == 8e3
        disp('Number of iterations too high')
        return
    else
%         disp('The relaxation method has achieved the required tolerance')
end
    
dLdr = 0.5*rho*Vef_nz.^2.*cl.*blade.chord;
dDdr = 0.5*rho*Vef_nz.^2.*cd.*blade.chord;
dTdr = (dLdr.*cos(beta-alpha_ef)-dDdr.*sin(beta-alpha_ef))./cos(blade.delta);
dQdr = (dLdr.*sin(beta-alpha_ef)+dDdr.*cos(beta-alpha_ef)).*r*blade.tip_rad;
T = blade.num*trapz(r*blade.tip_rad,dTdr);
P = omega*blade.num*trapz(r*blade.tip_rad,dQdr);
Ct = T/(rho*n^2*blade.diam^4);
Cp = P/(rho*n^3*blade.diam^5);
eta = T*V/P;
end