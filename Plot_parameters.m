function [] = Plot_parameters(mission,results,blade)
% Twist and chord
figure()
yyaxis left
plotting(results.r,blade.chord,'b','$\frac{r}{R}$','$c [m] $','Optimum chord and twist distribution')
hold on
yyaxis right
plotting(results.r,rad2deg(blade.beta),'r','$\frac{r}{R}$','$\beta [deg] $','Optimum chord and twist distribution')
legend('Chord','Twist')

% Performance Parameters
figure();
plotting_pars(results.adv_ratio,results.Ct,results.Cp,results.eta)
title('Performance parameters')
figure()
plotting(results.adv_ratio,results.P,'b','J','Power [W]','Power envelope')
hold on;plotting(0,results.P_static,'bo','J','Power [W]','Power envelope');
plotting(mission.j_cruise,results.P(end),'go','$J=\frac{V}{nd}$','Power [W]','Power envelope')
legend('Power','Static','Cruise')
figure();
plotting(results.adv_ratio,results.T,'b','J','Thrust [N]','Thrust envelope')
hold on;plotting(0,results.T_static,'bo','J','Thrust [N]','Thrust envelope');
plotting(mission.j_cruise,results.T(end),'go','$J=\frac{V}{nd}$','Thrust [N]','Thrust envelope')
legend('Thrust','Static','Cruise');hold off
% RPM's and pitch law
figure()
yyaxis left
plotting(results.adv_ratio,mission.omega*60/(2*pi),'b','$J=\frac{V}{nd}$','$\Omega [rpm]$','RPMs and VP laws')
yyaxis right
plotting(results.adv_ratio,rad2deg(results.pitch),'r','$J=\frac{V}{nd}$','$Pitch [deg]$','RPMs and VP laws')
legend('RPMs','Pitch')
function [] = plotting_pars(j,ct,cp,eta)
    yyaxis left
    grid on
    plot(j,ct,'LineWidth',1.5)
    hold on
    plot(j,cp,'LineWidth',1.5)
    xlabel('$J=\frac{V}{nd}$','Interpreter','latex');
    ylabel('$C_{t}$,$C_{p}$','Interpreter','latex');
    title('Performance parameters fixed rpms');
    yyaxis right
    plot(j,eta,'LineWidth',1.5)
    legend('$C_{t}$','$C_{p}$','$\eta$','Interpreter','latex')
    ylabel('$\eta$','Interpreter','latex')
end

function [] = plotting(x,y,col,xlab,ylab,tit)
plot(x,y,col,'LineWidth',1.5)
grid on
xlabel(xlab,'Interpreter','latex')
ylabel(ylab,'Interpreter','latex')
title(tit)
end

end