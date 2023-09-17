clear;
clc;
close all;

global I i j k
m = 5; % [kg] 
i = 2; % X distance from origin [m]
j = 1; % Y distance from origin [m]
k = 0.02; % Perturbation (Z) distance [m]

% Create matrix for mass/position for 4 point masses where each row
% corresponds to a point mass and the first 3 columns corresponds to x, y,
% z positions and the last column represents the mass.
mass_properties = [i 0 0 m; -i 0 0 m; 0 -j -k m/2; 0 j k m/2];
I = inertiaDyadic(mass_properties);

% Create vector to represent initial angular velocities p, q, r, phi,
% theta, psi
w0 = [2e-16 15 0 0 0 0]; % [rad/sec]

tstart = 0; % time start [sec]
tend = 30; % time end [sec]
n = 1000;
tspan = linspace(tstart, tend, n);

[t, angles] = ode23(@DzhanibekovIntFunction, tspan, w0);
angles = angles';

[q1, q2, q3, q4] = findPos(angles);

% Create graphs show p q r angles and Euler angles in one figure
f1 = figure(1);
f1.Position = [575 735 1200 500];
subplot(2,1,1)
plot(t, rad2deg(angles(1:3,:)))
title('Bx By Bz Angles')
legend('p', 'q', 'r')
qmax = rad2deg(max(angles(2,:)));
ylim([-1.25*qmax 1.25*qmax])

subplot(2,1,2)
plot(t, rad2deg(wrapTo2Pi(angles(4:6,:))))
title('Euler Angles')
legend('phi', 'theta', 'psi')
ylim([-10 375])

% Animation code, creates 3 different views for animation
f2 = figure(2);
f2.Position = [300 150 1750 500];
for h = 1:n
    for l = 1:3
    figure(2)
    subplot(1,3,l)
    plot3(q1(1,h), q1(2,h), q1(3,h), 'r.', 'markersize', 50)
    hold on
    plot3(q2(1,h), q2(2,h), q2(3,h), 'b.', 'markersize', 50)
    plot3(q3(1,h), q3(2,h), q3(3,h), 'g.', 'markersize', 20)
    plot3(q4(1,h), q4(2,h), q4(3,h), 'm.', 'markersize', 20)
    line([q1(1,h) q2(1,h)], [q1(2,h) q2(2,h)], [q1(3,h) q2(3,h)], 'Color', 'k', 'Linewidth', 3)
    line([q3(1,h) q4(1,h)], [q3(2,h) q4(2,h)], [q3(3,h) q4(3,h)], 'Color', 'k', 'Linewidth', 1)
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    xlim([-1.5*i 1.5*i])
    ylim([-1.5*i 1.5*i])
    zlim([-1.5*i 1.5*i])
    hold off
    grid on
    if l == 1
        view([0 0])
    elseif l == 2
        title('Dzhanibekov Effect Simulation')
    elseif l == 3
        view([-90 0])
    end
    end
end

function dy = DzhanibekovIntFunction(~, y)
    global I

    p = y(1);
    q = y(2);
    r = y(3);
    phi = y(4);
    theta = y(5);
    psi = y(6);
    w = [p; q; r];
    wSkew = [0 -r q;
        r 0 -p;
        -q p 0];
    dy(1:3) = -inv(I)*wSkew*I*w;
    dy(4) = p + (q*sin(phi) + r*cos(phi))*tan(theta);
    dy(5) = q*cos(phi) - r*sin(phi);
    dy(6) = (q*sin(phi) + r*cos(phi))*sec(theta);
    dy = dy';
end

function I = inertiaDyadic(m)
    % Creates Inertia Dyadic based on given parameters
    Ixx = m(3,4)*(m(3,2)^2 + m(3,3)^2) + m(4,4)*(m(4,2)^2 + m(4,3)^2);
    Iyy = m(1,4)*(m(1,1)^2 + m(1,3)^2) + m(2,4)*(m(2,1)^2 + m(2,3)^2) + m(3,4)*(m(3,1)^2 + m(3,3)^2) + m(4,4)*(m(4,1)^2 + m(4,3)^2);
    Izz = m(1,4)*(m(1,1)^2 + m(1,2)^2) + m(2,4)*(m(2,1)^2 + m(2,2)^2) + m(3,4)*(m(3,1)^2 + m(3,2)^2) + m(4,4)*(m(4,1)^2 + m(4,2)^2);
    Ixy = 0;
    Ixz = 0;
    Iyz = 0;
    I = [Ixx Ixy Ixz;
        Ixy Iyy Iyz;
        Ixz Iyz Izz];
end

function [Q1,Q2,Q3,Q4] = findPos(angles)
    global i j k

    phi = angles(4,:);
    theta = angles(5,:);
    psi = angles(6,:);

    x1 = i.*cos(theta).*cos(psi);
    y1 = i.*cos(theta).*sin(psi);
    z1 = -i.*sin(theta);

    x2 = -x1;
    y2 = -y1;
    z2 = -z1;

    x3 = -j*(sin(phi).*sin(theta).*cos(psi) - cos(phi).*sin(psi)) - k*(sin(phi).*sin(psi) + cos(phi).*sin(theta).*cos(psi));
    y3 = -j*(cos(phi).*cos(psi) + sin(phi).*sin(theta).*sin(psi)) - k*(cos(phi).*sin(theta).*sin(psi) - sin(phi).*cos(psi));
    z3 = -j.*sin(phi).*cos(theta) - k*cos(phi).*cos(theta);

    x4 = -x3;
    y4 = -y3;
    z4 = -z3;

    Q1 = [x1;y1;z1];
    Q2 = [x2;y2;z2];
    Q3 = [x3;y3;z3];
    Q4 = [x4;y4;z4];

end