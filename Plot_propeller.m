function [MP,LE,TE]=Plot_propeller(blade,results)

%Obtaining coordinates for LE and TE [x(radial),y(chordwise),z(twisted)]
% le = [results.r*blade.tip_rad.*cos(blade.sweep); blade.chord+results.r.*sin(blade.sweep)*blade.tip_rad; 0*results.r];
% te = [results.r*blade.tip_rad.*cos(blade.sweep); results.r.*sin(blade.sweep); 0*results.r];
% mp = [results.r*blade.tip_rad.*cos(blade.sweep); blade.chord/2+results.r.*sin(blade.sweep)*blade.tip_rad; 0*results.r];

le = [results.r*blade.tip_rad; blade.chord; 0*results.r];
te = [results.r*blade.tip_rad; results.r*0; 0*results.r];
mp = [results.r*blade.tip_rad; blade.chord/2; 0*results.r];

%Twist the blade obtaining the z coordinate
for i=1:blade.num
    angle = 2*pi/blade.num*i;
    rot_matrix_z = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];
    for j=1:length(results.r)
        rot_matrix_x = [1 0 0; 0 cos(blade.beta(j)) -sin(blade.beta(j)); 0 sin(blade.beta(j)) cos(blade.beta(j))];
        rot_matrix_sweep = [cos(-blade.sweep(j)) -sin(-blade.sweep(j)) 0; sin(-blade.sweep(j)) cos(-blade.sweep(j)) 0; 0 0 1];
        LE(:,j,i) = rot_matrix_z*(rot_matrix_sweep*(rot_matrix_x*le(:,j)));
        TE(:,j,i) = rot_matrix_z*(rot_matrix_sweep*(rot_matrix_x*te(:,j)));
        MP(:,j,i) = rot_matrix_z*(rot_matrix_sweep*(rot_matrix_x*mp(:,j)));
    end

end
% figure()
for i=1:blade.num
    plot3(LE(1,:,i),LE(2,:,i),LE(3,:,i),'k','LineWidth',2)
    hold on
    plot3(TE(1,:,i),TE(2,:,i),TE(3,:,i),'b','LineWidth',2)
    axis('equal')
    for j=1:length(results.r)
        plot3([LE(1,j,i) TE(1,j,i)],[LE(2,j,i) TE(2,j,i)],[LE(3,j,i) TE(3,j,i)],'k')
    end
end
