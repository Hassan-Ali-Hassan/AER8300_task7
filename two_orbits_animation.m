clear all
close all
clc


%% making a simple loop for calculating the circular trajectory

dt = 0.01;
omega_target = 1;
omega_orbit = 1;

r = 700;
alpha = 0; % alpha is the local angle in the local trajectory axes
beta = 0;
% R = 700000;
R = 1200;

% define the rotation matrix (constant)
rot_target = rotation(0,1,0);
rot_orbit = rotation(0, pi*0.5,0);

%the euler angles of the satellite itself
theta = 0;
phi = 0;
psi = 0;

for i = 1:10000 %which is the number of steps we want
    
    %calculating the local coorddinates wrt curcular trajectory plane
    x_local_target = [r*cos(alpha) , r*sin(alpha), 0]';
    x_local_orbit = [R*cos(beta) , R*sin(beta), 0]';
    
    %obtain the global coordinate of the point on trajectory
    x_global_target = rot_target  * x_local_target;
    x_global_orbit = rot_orbit  * x_local_orbit;
    
    %update the local angle alpha
    alpha = alpha + omega_target * dt;
    beta = beta + omega_orbit * dt;
    
    %save the global position
    X_target(i,:) = x_global_target';
    X_orbit(i,:) = x_global_orbit';
    
    %calculate the vector joining the target and the satellite
    x_rel_global = x_global_target - x_global_orbit;
    
    
end

%% For demonstration purpose
figure(1)
fh = plot3(X_target(:,1),X_target(:,2),X_target(:,3));
hold on
plot3(X_orbit(:,1),X_orbit(:,2),X_orbit(:,3));
axis equal

for i = 1:size(X_target,1)    
    fh = plot3(X_target(:,1),X_target(:,2),X_target(:,3));
    hold on
    plot3(X_orbit(:,1),X_orbit(:,2),X_orbit(:,3));
    plot3(X_target(i,1),X_target(i,2),X_target(i,3),'bo','MarkerSize',8,'MarkerFaceColor','b')
    plot3(X_orbit(i,1),X_orbit(i,2),X_orbit(i,3),'ro','MarkerSize',8,'MarkerFaceColor','r')
    line([X_orbit(i,1),X_target(i,1)],[X_orbit(i,2),X_target(i,2)],[X_orbit(i,3),X_target(i,3)],'Color','green')
    drawnow
    cla
end



%% 
function Rt = rotation(phi,theta,psi)
Rx = [1 0 0;0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta) ];
Rz = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];

Rt = Rx * Ry* Rz; %target rotation matrix
end

