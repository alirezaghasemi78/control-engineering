%% ACTUATOR SATURATION

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

x0 = [30*pi/180;0];

u_limits = [Inf 5 2 1 0.5];

figure;

for k = 1:length(u_limits)

    umax = u_limits(k);

    x = zeros(2,length(t));
    u = zeros(1,length(t));

    x(:,1) = x0;

    for i = 1:length(t)-1

        u_unsat = -K*x(:,i);

        if isfinite(umax)

            u(i) = min(max(u_unsat,-umax),umax);

        else

            u(i) = u_unsat;

        end

        xdot = A*x(:,i)+B*u(i);

        x(:,i+1) = x(:,i)+dt*xdot;

    end

    plot(t,x(1,:),'LineWidth',1.3);
    hold on;

end

grid on;

xlabel('Time [s]');
ylabel('\theta [rad]');

legend( ...
    'No saturation',...
    '|u|<=5',...
    '|u|<=2',...
    '|u|<=1',...
    '|u|<=0.5');

title('Effect of Actuator Saturation');