%% ============================================================
% NONLINEAR PLANT vs LINEAR PLANT
% LQR designed from linearized model
% ============================================================

clear;
clc;
close all;

%% Parameters

m = 1;
l = 1;
b = 0.2;
g = 9.81;

%% ============================================================
% 1. Linearized model
% ============================================================

A = [0 1;
     g/l -b/(m*l^2)];

B = [0;
     1/(m*l^2)];

C = [1 0];
D = 0;

%% ============================================================
% 2. Check open-loop stability
% ============================================================

eig_A = eig(A);

disp('==========================================');
disp('Open-loop linear eigenvalues:');
disp(eig_A);

%% ============================================================
% 3. LQR controller
% ============================================================

Q = diag([100 10]);
R = 1;

K = lqr(A,B,Q,R);

disp(' ');
disp('LQR gain K:');
disp(K);

%% Closed-loop linear system

Acl = A - B*K;

eig_Acl = eig(Acl);

disp(' ');
disp('Closed-loop linear eigenvalues:');
disp(eig_Acl);

%% ============================================================
% 4. Initial conditions
% ============================================================

initial_angles_deg = [5 15 45 60 75 90];

tspan = [0 5];

%% Storage

results = struct();

%% ============================================================
% 5. Simulation
% ============================================================

for k = 1:length(initial_angles_deg)

    theta0 = initial_angles_deg(k)*pi/180;

    x0 = [theta0; 0];

    %% --------------------------------------------------------
    % Linear plant
    % --------------------------------------------------------

    [t_lin,x_lin] = ode45( ...
        @(t,x) linear_pendulum(t,x,K,A,B), ...
        tspan, ...
        x0);

    %% --------------------------------------------------------
    % Nonlinear plant
    % --------------------------------------------------------

    [t_non,x_non] = ode45( ...
        @(t,x) nonlinear_pendulum(t,x,K,m,l,b,g), ...
        tspan, ...
        x0);

    %% --------------------------------------------------------
    % Control signals
    % --------------------------------------------------------

    u_lin = zeros(length(t_lin),1);
    u_non = zeros(length(t_non),1);

    for i = 1:length(t_lin)
        u_lin(i) = -K*x_lin(i,:)';
    end

    for i = 1:length(t_non)
        u_non(i) = -K*x_non(i,:)';
    end

    %% Store results

    results(k).angle = initial_angles_deg(k);

    results(k).t_lin = t_lin;
    results(k).x_lin = x_lin;
    results(k).u_lin = u_lin;

    results(k).t_non = t_non;
    results(k).x_non = x_non;
    results(k).u_non = u_non;

end

%% ============================================================
% 6. Plot 1
% Linear vs Nonlinear response
% ============================================================

figure;

for k = 1:length(initial_angles_deg)

    subplot(3,2,k);

    plot(results(k).t_lin, ...
         results(k).x_lin(:,1)*180/pi, ...
         'LineWidth',1.5);

    hold on;

    plot(results(k).t_non, ...
         results(k).x_non(:,1)*180/pi, ...
         '--', ...
         'LineWidth',1.5);

    grid on;

    xlabel('Time [s]');
    ylabel('\theta [deg]');

    title(sprintf('Initial angle = %d deg', ...
        initial_angles_deg(k)));

    legend('Linear plant','Nonlinear plant');

end

sgtitle('Linear vs Nonlinear Plant with Same LQR Controller');

%% ============================================================
% 7. Plot 2
% Difference between linear and nonlinear models
% ============================================================

figure;

for k = 1:length(initial_angles_deg)

    theta_lin = interp1( ...
        results(k).t_lin, ...
        results(k).x_lin(:,1), ...
        results(k).t_non);

    error = results(k).x_non(:,1) - theta_lin;

    subplot(3,2,k);

    plot(results(k).t_non, ...
         error*180/pi, ...
         'LineWidth',1.5);

    grid on;

    xlabel('Time [s]');
    ylabel('\Delta\theta [deg]');

    title(sprintf('Initial angle = %d deg', ...
        initial_angles_deg(k)));

end

sgtitle('Difference Between Nonlinear and Linear Responses');

%% ============================================================
% 8. Plot 3
% Control effort
% ============================================================

figure;

for k = 1:length(initial_angles_deg)

    subplot(3,2,k);

    plot(results(k).t_lin, ...
         results(k).u_lin, ...
         'LineWidth',1.5);

    hold on;

    plot(results(k).t_non, ...
         results(k).u_non, ...
         '--', ...
         'LineWidth',1.5);

    grid on;

    xlabel('Time [s]');
    ylabel('u');

    title(sprintf('Initial angle = %d deg', ...
        initial_angles_deg(k)));

    legend('Linear','Nonlinear');

end

sgtitle('Control Effort');

%% ============================================================
% 9. Quantitative comparison
% ============================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('        LINEAR vs NONLINEAR COMPARISON\n');
fprintf('====================================================\n');

fprintf('%10s %15s %15s %15s\n', ...
    'Angle','Max Error','Final NL','Max Control');

fprintf('----------------------------------------------------\n');

for k = 1:length(initial_angles_deg)

    theta_lin = interp1( ...
        results(k).t_lin, ...
        results(k).x_lin(:,1), ...
        results(k).t_non);

    theta_non = results(k).x_non(:,1);

    error = theta_non - theta_lin;

    max_error = max(abs(error))*180/pi;

    final_theta = theta_non(end)*180/pi;

    max_control = max(abs(results(k).u_non));

    fprintf('%10d %15.4f %15.4f %15.4f\n', ...
        initial_angles_deg(k), ...
        max_error, ...
        final_theta, ...
        max_control);

end

fprintf('====================================================\n');


%% ============================================================
% LINEAR PLANT
% ============================================================

function dx = linear_pendulum(~,x,K,A,B)

    u = -K*x;

    dx = A*x + B*u;

end


%% ============================================================
% NONLINEAR PLANT
% ============================================================

function dx = nonlinear_pendulum(~,x,K,m,l,b,g)

    theta = x(1);
    theta_dot = x(2);

    %% LQR controller

    u = -K*x;

    %% Nonlinear dynamics

    theta_ddot = ...
        (g/l)*sin(theta) ...
        - (b/(m*l^2))*theta_dot ...
        + (1/(m*l^2))*u;

    dx = [theta_dot;
          theta_ddot];

end