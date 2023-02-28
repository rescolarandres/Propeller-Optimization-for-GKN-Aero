%%--BLADE OPTIMIZATION CODE USING BEMT & VORTEX THEORY BY RODRIGO ESCOLAR ANDRES/November 2021
clc;clear;close all
% 
%% Input data
blade.diam = 2.8; %[m]     BLADE DIAMETER
blade.tip_rad = blade.diam/2; %[m]   TIP RADIUS
blade.hub_rad = blade.tip_rad*0.2;%[m]   HUB RADIUS
blade.num = 5; %[num_blades]     NUMBER OF BLADES
blade.alpha = deg2rad(8); %[rad]     DESIRED AoA AT BLADE ELEMENT (obtain min Cd/Cl ratio for airfoil)
blade.cd = importdata('CDCLARKY.csv');    % AERO COEFFS FOR CHOSEN AIRFOIL
blade.cl = importdata('CLCLARKY.csv');
blade.radial_sections = 15;
blade.sweep = linspace(0,deg2rad(15),blade.radial_sections);

load('t_req_climb.mat');load('optim_results_climb.mat');
mission.optim = optim_results;
mission.flight_step = 15;   % NUMBER OF DIVISIONS OF FLIGHT ENVELOPE
mission.mach = 0.5;     % CRUISE MACH NUMBER
mission.omega = linspace(2100,1100,mission.flight_step)*2*pi/60; %[rad/s]   NON OPTIMIZED RPM'S LAW FOR WHOLE FLIGHT
mission.T_req = 4000;  %[N]    THRUST REQUIRED AT OPTIMIZATION CONDITION
mission.t_req =  t_req; %THRUST REQUIRED DURING CLIMBING
mission.altitude = linspace(0,25000,mission.flight_step)*0.3048; %[m]    ALTITUDE ENVELOPE
mission.viscosity = 1.789e-5;%  VISCOSITY FLUID (AIR)
[mission.T,mission.a,mission.P,mission.rho] = atmosisa(mission.altitude); %[SI]   THERMO PARAMETERS FOR FLIGHT ENVELOPE
mission.V = mission.mach*mission.a(mission.flight_step); %[m/s]     CRUISE VELOCITY
mission.V_takeoff = 66.88;%[m/s]    TAKE OFF VELOCITY
mission.j_cruise = mission.V/(mission.omega(mission.flight_step)*blade.diam/(2*pi));%   CRUISE ADVANCE RATIO
mission.j_takeoff = mission.V_takeoff/(mission.omega(1)/(2*pi)*blade.diam);%   TAKE OFF ADVANCE RATIO
mission.Re = mission.rho(mission.flight_step)*mission.V*0.1/mission.viscosity;%     REYNOLDS AT OPTIMIZATION CONDITION

%% Optimization process to obtain blade geometry
results.r = linspace(blade.hub_rad,blade.tip_rad,length(blade.sweep))/blade.tip_rad; %   NON DIM RADIUS
results.pitch = linspace(deg2rad(-6),deg2rad(40),mission.flight_step);%  NON OPTIMIZED VAR PITCH LAW
results.adv_ratio = linspace(0,mission.j_cruise,mission.flight_step);%  FLIGHT ENVELOPE ADV RATIO
blade.delta = pi/2-atan(diff(results.r*blade.tip_rad.*cos(blade.sweep))./(diff(results.r*blade.tip_rad.*sin(blade.sweep))));    % Blade sweep derivative distribution
blade.delta = [0 blade.delta];

x = 2*linspace(blade.hub_rad,blade.tip_rad,length(results.r))/blade.diam;% NON DIM X
% Blade geometry is obtained using solve optimum function, by specifying design conditions: omega, V and rho as a function of design variables: AoA, D, N
[blade.chord,blade.beta,results.a,results.b] = solve_optimum (mission.omega(mission.flight_step),mission.rho(mission.flight_step),mission,blade,results.r,x,'v2');
% Removing the pitch from blade twist as the optimizer computes total beta
blade.beta = blade.beta-results.pitch(length(results.pitch));
% Plot optimized propeller
[blade.MP,blade.LE,blade.TE]=Plot_propeller(blade,results);

%% Obtaining performances as a function of advance ratio J for a given blade
% As the blade geometry is given now, prop is analyzed using BEMT as a function of the design variables

f_bemt = @(j,beta,omega,rho) BEMT_RE(j,mission,blade,results.r,beta,omega,rho);
f_vortex = @(j,beta,omega,rho) VORTEX_RE(j,mission,blade,results.r,beta,omega,rho);

% % Take off performace
% load('takeoff.mat')
% for i=1:length(takeoff.velocities)
%     takeoff.pitch(i)= optimize_takeoff(takeoff.velocities(i)/(blade.diam*2100/60),mission,blade,results,mission.rho(1),2100*0.1047);
%     [~,~,~,takeoff.T(i),takeoff.P(i),takeoff.Ve(i,:),~] = f_vortex(takeoff.velocities(i)/(blade.diam*2100/60),(blade.beta+takeoff.pitch(i)),2100*0.1047,mission.rho(1));
% 
% end

% Climbing performance
% for i=1:length(results.adv_ratio)
% % %     optim_results(i,:) = optimize_climb(results.adv_ratio(i),mission,blade,results,mission.rho(i),i);
%     results.pitch(i)=optim_results(i,1);
%     mission.omega(i)=optim_results(i,2);
%     [results.Ct(i),results.Cp(i),results.eta(i),results.T(i),results.P(i),results.Ve(i,:),~,~] = f_vortex(results.adv_ratio(i),(blade.beta+results.pitch(i)),mission.omega(i),mission.rho(i));
% end

% 
% %  Performance for various conditions
% [results.Ct_static,results.Cp_static,results.eta_static,results.T_static,results.P_static,results.V,~] = f_vortex(results.adv_ratio(1),blade.beta+results.pitch(1),mission.omega(1),mission.rho(1));
% [results.Ct_cruise,results.Cp_cruise,results.eta_cruise,results.T_cruise,results.P_cruise,~,~,~,results.w_cruise] = f_vortex(mission.j_cruise,(blade.beta+results.pitch(mission.flight_step)),mission.omega(mission.flight_step),mission.rho(mission.flight_step));
% results.mach = (results.V(end))./mission.a(1);

% % Descent performance
% load('descent_700_low.mat')
% for i=1:length(descent.altitudes)
%     f=@(pitch) optimize_descent(descent.vel/(blade.diam*descent.omega*2*pi),mission,blade,results.r,pitch,descent.rho(i),descent.omega,mission.t_req_descent(i));
%     descent.pitch(i)=fzero(f,0.3,optimset('Display','iter'));
%     [descent.Ct(i),descent.Cp(i),descent.eta(i),descent.T(i),descent.P(i),~,~] = f_vortex(descent.j(i),(blade.beta+descent.pitch(i)),descent.omega,descent.rho(i));
% end
% 
figure()
yyaxis right
plot(results.adv_ratio,rad2deg(optim_results(:,1)),'r','LineWidth',1.5);hold on
ylabel('Pitch [deg]')
yyaxis left
plot(results.adv_ratio,optim_results(:,2)*60/(2*pi),'b','LineWidth',1.5);hold on
ylabel('$\Omega$ [rpm]','Interpreter','latex')
xlabel('$J=\frac{V}{nd}$','Interpreter','latex')
legend('Variable pitch setting','Rotational speed setting')
title('Climb settings')
grid on
%% Plotting parameters

Plot_parameters(mission,results,blade);

function [chord,twist,a,b] = solve_optimum (omega,rho,mission,blade,r,x,version)
    % This function solves blade geometry by means of a NLS until the desired thrust is obtained. There is currently two optimization methods
    if version == 'v1'
        % This method solves induction axial velocity by iterating BEMT
        f = @(v0) shootv1(v0,mission,blade,r,x,omega,rho);
        v0 = fzero(f,10);
        % As the value of v0 that gives that thrust for minimal power is obtained, the code is run to obtain the final geometry
        [~,chord,twist,a,b] = shootv1 (v0,mission,blade,r,x,omega,rho);
    
    elseif version == 'v2'
        % This method solves the vortex tip velocity by iterating BEMT
        f = @(chi) shootv2(chi,mission,blade,r,omega,rho);
        chi = fzero(f,0);
        [~,chord,twist,a,b] = shootv2 (chi,mission,blade,r,omega,rho);
    end
    function [error,chord,beta,alpha,a,b] = shootv1 (v0,mission,blade,r,x,omega,rho)
    phi = atan((mission.V+v0)./(omega*r));
    [cl,cd] = get_coeffs(blade.alpha*180/pi,mission.Re,blade);
    epsilon = cd/cl;
    a = v0*cos(phi).^2/mission.V.*(1-epsilon*tan(phi));
    b = v0./(omega*r).*cos(phi).*sin(phi).*(1+epsilon./tan(phi));
    Vr = sqrt((mission.V*(1+a)).^2+(omega*r.*(1-b)).^2);
    %--Tip loss factor calculation
    F = 2/pi.*acos(exp(blade.num*omega*blade.diam.*(x-1)/(4*mission.V)));
    beta = blade.alpha + phi;
    %--Circulation is computed so the solidity can be obtained
    circ = pi*x.^2.*b*omega*blade.diam^2.*F/blade.num;
    sol = 2/pi*blade.num*circ./(cl*Vr.*x*blade.diam);
    chord = 2*pi*r.*sol/blade.num;         
                 
    [~,~,~,T,~,~,~,~,alpha]=BEMT_RE(mission.j_cruise,mission,blade,r,beta,chord,omega,rho);
    % Check computed T with required one, if it is less than tolerance iterate again
    error = T-mission.T_req;
    end
function [error,chord,beta,a,b] = shootv2 (chi,mission,blade,r,omega,rho)

    lambda = mission.V/(omega*blade.tip_rad);
    phi_t = atan(lambda*(1+chi/2));
    phi = atan(tan(phi_t)./r);
    F = 2/pi*acos(exp(-blade.num*0.5*(1-r)/sin(phi_t)));
    G = F.*cos(phi).*sin(phi);
    [cl,cd] = get_coeffs(rad2deg(blade.alpha),mission.Re,blade);
    Wc = 4*pi*lambda*G*mission.V*blade.tip_rad*chi/(cl*blade.num);
    epsilon = cd/cl;
    a = chi*0.5*cos(phi).^2.*(1-epsilon*tan(phi));
    b = chi./(2*omega.*r/mission.V).*cos(phi).*sin(phi).*(1+epsilon./tan(phi));
    W = mission.V*(1+a)./sin(phi);
    chord = Wc./W;
    beta = blade.alpha+phi;
    blade.chord = chord;

    [~,~,~,T,~]=BEMT_RE(mission.j_cruise,mission,blade,r,beta,omega,rho);
    error = T-mission.T_req;
end
    end

function [x] = optimize_climb(j,mission,blade,results,rho,i)
    % This function obtains the optimum pitch as a function of max/min
    % power or thrust 
    f_pitch = @(x) opt_pitch(j,mission,blade,results.r,x,rho);
    nonlcon = @(x) constraint(j,mission,blade,results.r,x,rho,i);
    if i ~= 1
        x0 = [mission.optim(i-1,1), mission.optim(i-1,2)];
        lb = [(results.pitch(i-1)-deg2rad(15)), (mission.omega(i-1)-500*0.1047)];
        ub = [(results.pitch(i-1)+deg2rad(15)), (mission.omega(i-1)+100*0.1046)];
    else
        x0 = [results.pitch(1), mission.omega(1)];
        lb = [(results.pitch(1)-deg2rad(10)), (mission.omega(1)-500*0.1047)];
        ub = [(results.pitch(1)+deg2rad(10)), (mission.omega(1)+100*0.1046)];
    end
    x = fmincon(f_pitch,x0,[],[],[],[],lb,ub,nonlcon,optimoptions('fmincon','Display','iter'));

    function [P] = opt_pitch(j,mission,blade,r,x,rho)
       % This function obtains T or P depending on what we want to max/min
    [~,~,~,~,P,~,~,~,~]=VORTEX_RE(j,mission,blade,r,x(1)+blade.beta,x(2),rho);
    end

    function [c,ceq] = constraint(j,mission,blade,r,x,rho,i)
        % This function gives the constrains in power usage or thrust
        % required
        [~,~,~,T,~,~,~,~,~]=VORTEX_RE(j,mission,blade,r,x(1)+blade.beta,x(2),rho);
        c = mission.t_req_descent(i)-T;
        ceq = [];
    end
end

function[error] = optimize_descent(j,mission,blade,r,pitch,rho,omega,t_req)
    [~,~,~,T,~,~,~,~,~]=VORTEX_RE(j,mission,blade,r,pitch+blade.beta,omega,rho);
    error = T-t_req;
end

function [x] = optimize_takeoff(j,mission,blade,results,rho,omega)
% This function obtains the optimum pitch as a function of max/min
    % power or thrust 
    f_pitch = @(pitch) -opt_pitch(j,mission,blade,results.r,omega,pitch,rho);
    nonlcon = @(pitch) constraint(j,mission,blade,results.r,omega,pitch,rho);
    x = fmincon(f_pitch,-0.0759,[],[],[],[],[deg2rad(-10)],[deg2rad(10)],nonlcon,optimoptions('fmincon','Display','iter'));

    function [T] = opt_pitch(j,mission,blade,r,omega,pitch,rho)
       % This function obtains T or P depending on what we want to max/min
    [~,~,~,T,~,~,~,~,~]=VORTEX_RE(j,mission,blade,r,pitch+blade.beta,omega,rho);
    end
    
    function [c,ceq] = constraint(j,mission,blade,r,omega,pitch,rho)
        % This function gives the constrains in power usage or thrust
        % required
        [~,~,~,~,P,~,~,~,~]=VORTEX_RE(j,mission,blade,r,pitch+blade.beta,omega,rho);
        c = P-1e6;
        ceq = [];
    end


end
