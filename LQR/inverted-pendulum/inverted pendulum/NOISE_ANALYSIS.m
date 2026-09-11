%% NOISE ANALYSIS

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
t = 0:dt:10;

x0 = [10*pi/180;0];

noise_levels = [0 0.001 0.005 0.01 0.05 0.1];

figure;

for k = 1:length(noise_levels)

    sigma = noise_levels(k);

    x = zeros(2,length(t));
    u = zeros(1,length(t));

    x(:,1) = x0;

    for i = 1:length(t)-1

        noise = sigma*randn(2,1);

        x_measured = x(:,i) + noise;

        u(i) = -K*x_measured;

        xdot = A*x(:,i) + B*u(i);

        x(:,i+1) = x(:,i) + dt*xdot;

    end

    plot(t,x(1,:),'LineWidth',1.2);
    hold on;

end

grid on;

xlabel('Time [s]');
ylabel('\theta [rad]');

legend( ...
    '\sigma=0',...
    '\sigma=0.001',...
    '\sigma=0.005',...
    '\sigma=0.01',...
    '\sigma=0.05',...
    '\sigma=0.1');

title('Effect of Measurement Noise');