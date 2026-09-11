%% DISTURBANCE REJECTION

clear;
clc;
close all;

%% Plant

m = 1;
l = 1;
b = 0.2;
g = 9.81;

A = [0 1;
     g/l -b/(m*l^2)];

B = [0;
     1/(m*l^2)];

%% LQR

Q = diag([100 10]);
R = 1;

K = lqr(A,B,Q,R);

%% Simulation

dt = 0.001;
t = 0:dt:5;

x0 = [5*pi/180;0];

disturbance_levels = [0 0.1 0.5 1];

figure;

for k = 1:length(disturbance_levels)

    d_amp = disturbance_levels(k);

    x = zeros(2,length(t));

    x(:,1) = x0;

    for i = 1:length(t)-1

        % External disturbance

        if t(i) >= 1 && t(i) <= 2
            d = d_amp;
        else
            d = 0;
        end

        u = -K*x(:,i);

        xdot = A*x(:,i) + B*(u+d);

        x(:,i+1) = x(:,i) + dt*xdot;

    end

    plot(t,x(1,:),'LineWidth',1.3);
    hold on;

end

grid on;

xlabel('Time [s]');
ylabel('\theta [rad]');

legend( ...
    'd=0',...
    'd=0.1',...
    'd=0.5',...
    'd=1');

title('LQR Disturbance Rejection');