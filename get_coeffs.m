function [Cl,Cd]=get_coeffs(alpha,Re,blade)
% This function obtains the aerodynamic coefficients from the loaded data
% Given the Reynolds number and the angle of attack, the coefficients are
% obtained by a slpline interpolation 
% Input alpha must be in degrees
alphas = blade.cl(:,1);
cls = blade.cl(:,2:11);     %0.64 with ice
cds = blade.cd(:,2:11);     %1.25 with ice
ReInt = [60000 80000 100000 130000 160000 200000 300000 500000 1000000 3000000];

% 1. Get the Reynolds closest to the set of Reynolds available (to avoid 2D interpolation) as Reynolds effects are of minor importance
% 2. Interpolate on AoA to obtain coefficient
[~, idx] = min(abs(ReInt-Re));
Cl = spline(alphas,cls(:,idx),alpha);
Cd = spline(alphas,cds(:,idx),alpha);

end