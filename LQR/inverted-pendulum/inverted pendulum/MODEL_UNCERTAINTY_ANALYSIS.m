%% MODEL UNCERTAINTY ANALYSIS

clear;
clc;
close all;

%% Nominal model

m_nom = 1;
l_nom = 1;
b_nom = 0.2;
g = 9.81;

A_nom = [0 1;
         g/l_nom -b_nom/(m_nom*l_nom^2)];

B_nom = [0;
         1/(m_nom*l_nom^2)];

%% LQR designed for nominal model

Q = diag([100 10]);
R = 1;

K = lqr(A_nom,B_nom,Q,R);

%% Actual plant variations

mass_values = [0.5 0.75 1 1.25 1.5];

x0 = [10*pi/180;0];

dt = 0.001;
t = 0:dt:5;

figure;

for k = 1:length(mass_values)

    m = mass_values(k);

    A = [0 1;
         g/l_nom -b_nom/(m*l_nom^2)];

    B = [0;
         1/(m*l_nom^2)];

    x = zeros(2,length(t));

    x(:,1) = x0;

    for i = 1:length(t)-1

        u = -K*x(:,i);

        xdot = A*x(:,i)+B*u;

        x(:,i+1) = x(:,i)+dt*xdot;

    end

    plot(t,x(1,:),'LineWidth',1.3);
    hold on;

end

grid on;

xlabel('Time [s]');
ylabel('\theta [rad]');

legend( ...
    'm=0.5',...
    'm=0.75',...
    'm=1',...
    'm=1.25',...
    'm=1.5');

title('LQR Robustness to Mass Uncertainty');