function [V_ef_etaz,alpha_ef,cl,cd,w,alpha_i,circ] = circulation_sweep(alpha_i0,circ0,V,r,blade,mission,beta,omega,rho)
% This function obtains the circulation along the blades numercially,
% considering swept blades according to Tremmel's paper
r = r*blade.tip_rad; 
dphi = deg2rad(1);    % [rad] Discretization of vortex filament
% 1. Computation of w substituting integrals for infenitesimal sums
phi = atan(V./(omega*r));
w_cart = zeros(length(r),3);

for h=1:length(r)
    wp = 0; % Bound vortex integral
    wt = 0; % Trailing vortex integral
    ref = blade.MP(:,h,1)';

    for j = 1:blade.num
        % Induced velocity by bound vortices
        dwp = zeros(length(r)-1,3);
        dwt = zeros(length(r)-1,3);
        for i=2:length(r)
            dlp = blade.MP(:,i,j)'-blade.MP(:,i-1,j)';
            lp = ((blade.MP(:,i,j)+blade.MP(:,i-1,j))*0.5)';
            sp = ref-lp;
            dwp(i,:) = 0.5*(circ0(i)+circ0(i-1))*cross(dlp,sp)./(norm(sp)^3*4*pi);
        
%         Induced velocity by trailing vortices
            lt_p(1,:) = [(r(i)+r(i-1))/2 2*pi/blade.num*j-dphi/2 0]; % Polar coordinates
            lt_c(1,:) = lp;
            int = 0;
            for k=2:360
                dlt_p = [0 -dphi -sin(phi(j)+alpha_i0(j))*r(j)*dphi];
                lt_p(k,:) = lt_p(k-1,:) + dlt_p;
                lt_c(k,:) = [lt_p(k,1)*cos(lt_p(k,2)) lt_p(k,1)*sin(lt_p(k,2)) lt_p(k,3)];
                dlt_c = lt_c(k,:)-lt_c(k-1,:);
                st = ref-lt_c(k,:);
                int = int + cross(dlt_c,st)/norm(st)^3;
            end
            dgamma = 0.5*(circ0(i)+circ0(i-1));
            dwt(i,:) = -dgamma/(4*pi)*int;
%             plot3(lt_c(:,1),lt_c(:,2),lt_c(:,3),'k')
%             hold on
        end
        dwp(isnan(dwp))=0;
        wp =  wp + trapz(dwp);
        wt = wt + trapz(r(1)-r(2),dwt);
    end
    w_cart(h,:) = wp+wt;
end

clearvars -except w_cart r alpha_i0 circ0 phi blade mission V omega rho beta
[theta,~] = cart2pol(blade.MP(1,:,1),blade.MP(2,:,1),blade.MP(3,:,1));

w(:,1) = w_cart(:,1).*cos(theta')+w_cart(:,2).*sin(theta');
w(:,2) = -w_cart(:,1).*sin(theta')+w_cart(:,2).*cos(theta');
w(:,3) = w_cart(:,3);

% 2. Computing resultant velocity, alpha, Re for each radial section
alpha_i = 0.*r; V_ef_eta = 0.*r; V_ef_etaz = 0.*r; cl = 0.*r; cd = cl; circ = cl; alpha_ef = alpha_i;

for l=1:length(r)
    V_ef(l,:) = w(l,:)+[0 -omega*r(l) -V];
    alpha_i(l) = atan(-V_ef(l,3)/-V_ef(l,2))-phi(l);
    V_ef_eta(l) = -sin(blade.delta(l))*V_ef(l,1)+cos(blade.delta(l))*V_ef(l,2);
    alpha_ef(l) = beta(l)-atan(-V_ef(l,3)/-V_ef_eta(l));
    V_ef_etaz(l) = sqrt(V_ef_eta(l)^2+V_ef(l,3)^2);
    Re = V_ef_etaz(l)*blade.chord(l)/(mission.viscosity/rho);
    [cl(l),cd(l)]=get_coeffs(rad2deg(alpha_ef(l)),Re,blade);
    circ(l) = 0.5*V_ef_etaz(l)*cl(l)*blade.chord(l);
end

% 3. Boundary conditions are enforced
circ(length(r))=0; circ(1)=0;

end


