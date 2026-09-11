%% ============================================================
% NONLINEAR CART-POLE
% LQR designed from LINEARIZED model
% Controller applied to NONLINEAR plant
% ============================================================

clear;
clc;
close all;

%% ============================================================
% 1. Parameters
% ============================================================

M = 1.0;       % Cart mass [kg]
m = 0.2;       % Pendulum mass [kg]
l = 0.5;       % Distance to pendulum center of mass [m]
b = 0.1;       % Cart friction
g = 9.81;      % Gravity

%% ============================================================
% 2. Linearized model around upright equilibrium
% ============================================================

A = [0                  1              0                         0;
     0                 -b/M         m*g/M                     0;
     0                  0              0                         1;
     0              -b/(M*l)       (M+m)*g/(M*l)                0];

B = [0;
     1/M;
     0;
    1/(M*l)];

%% ============================================================
% 3. Check open-loop stability
% ============================================================

disp('==============================================');
disp('OPEN LOOP EIGENVALUES');
disp('==============================================');

eig_open = eig(A)

%% ============================================================
% 4. LQR Controller
% ============================================================

Q = diag([10 1 100 10]);
R = 0.1;

K = lqr(A,B,Q,R);

disp('==============================================');
disp('LQR GAIN');
disp('==============================================');

K

%% ============================================================
% 5. Closed-loop linear stability
% ============================================================

Acl = A - B*K;

disp('==============================================');
disp('CLOSED LOOP LINEAR EIGENVALUES');
disp('==============================================');

eig_closed = eig(Acl)

%% ============================================================
% 6. Initial angles
% ============================================================

initial_angles_deg = [2 10];

tspan = [0 5];

%% ============================================================
% 7. Simulation
% ============================================================

results = struct();

for k = 1:length(initial_angles_deg)

    theta0 = initial_angles_deg(k)*pi/180;

    % Initial state:
    %
    % x = [position;
    %      velocity;
    %      angle;
    %      angular velocity]

    x0 = [0;
          0;
          theta0;
          0];

    %% --------------------------------------------------------
    % LINEAR PLANT
    % --------------------------------------------------------

    [t_lin,x_lin] = ode45( ...
        @(t,x) linear_cartpole(t,x,A,B,K), ...
        tspan, ...
        x0);

    %% --------------------------------------------------------
    % NONLINEAR PLANT
    % --------------------------------------------------------

    [t_non,x_non] = ode45( ...
        @(t,x) nonlinear_cartpole(t,x,K,M,m,l,b,g), ...
        tspan, ...
        x0);

    %% --------------------------------------------------------
    % Control signals
    % --------------------------------------------------------

    % Linear control input
    u_lin = -(K*x_lin')';

    % Nonlinear control input
    u_non = -(K*x_non')';

    %% Store

    results(k).angle = initial_angles_deg(k);

    results(k).t_lin = t_lin;
    results(k).x_lin = x_lin;
    results(k).u_lin = u_lin;

    results(k).t_non = t_non;
    results(k).x_non = x_non;
    results(k).u_non = u_non;

end

%% ============================================================
% Plot settings
% ============================================================

font_size = 12;
title_size = 13;
line_width = 1.8;

%% ============================================================
% FIGURE 1
% Pendulum Angle
% ============================================================

fig1 = figure( ...
    'Units','normalized', ...
    'Position',[0.05 0.05 0.90 0.85], ...
    'Color','w');

for k = 1:length(initial_angles_deg)

    subplot(1,length(initial_angles_deg),k);

    plot(results(k).t_lin, ...
         results(k).x_lin(:,3)*180/pi, ...
         'LineWidth',line_width);

    hold on;

    plot(results(k).t_non, ...
         results(k).x_non(:,3)*180/pi, ...
         '--', ...
         'LineWidth',line_width);

    yline(0,'k:','LineWidth',1.2);

    grid on;
    box on;

    xlabel('Time [s]', ...
        'FontSize',font_size);

    ylabel('\theta [deg]', ...
        'FontSize',font_size);

    title(sprintf('Initial angle = %d deg', ...
        initial_angles_deg(k)), ...
        'FontSize',title_size);

    legend('Linear','Nonlinear', ...
        'Location','best', ...
        'FontSize',10);

    set(gca, ...
        'FontSize',font_size, ...
        'LineWidth',1.0);

end

sgtitle('Linear vs Nonlinear Cart-Pole: Pendulum Angle', ...
    'FontSize',15, ...
    'FontWeight','bold');

%% ============================================================
% FIGURE 2
% Cart Position
% ============================================================

fig2 = figure( ...
    'Units','normalized', ...
    'Position',[0.05 0.05 0.90 0.85], ...
    'Color','w');

for k = 1:length(initial_angles_deg)

    subplot(1,length(initial_angles_deg),k);

    plot(results(k).t_lin, ...
         results(k).x_lin(:,1), ...
         'LineWidth',line_width);

    hold on;

    plot(results(k).t_non, ...
         results(k).x_non(:,1), ...
         '--', ...
         'LineWidth',line_width);

    grid on;
    box on;

    xlabel('Time [s]', ...
        'FontSize',font_size);

    ylabel('Cart Position [m]', ...
        'FontSize',font_size);

    title(sprintf('Initial angle = %d deg', ...
        initial_angles_deg(k)), ...
        'FontSize',title_size);

    legend('Linear','Nonlinear', ...
        'Location','best', ...
        'FontSize',10);

    set(gca, ...
        'FontSize',font_size, ...
        'LineWidth',1.0);

end

sgtitle('Linear vs Nonlinear Cart-Pole: Cart Position', ...
    'FontSize',15, ...
    'FontWeight','bold');

%% ============================================================
% FIGURE 3
% Control Effort
% Linear vs Nonlinear
% ============================================================

fig3 = figure( ...
    'Units','normalized', ...
    'Position',[0.05 0.05 0.90 0.85], ...
    'Color','w');

for k = 1:length(initial_angles_deg)

    subplot(1,length(initial_angles_deg),k);

    % Linear control signal
    plot(results(k).t_lin, ...
         results(k).u_lin, ...
         'LineWidth',line_width);

    hold on;

    % Nonlinear control signal
    plot(results(k).t_non, ...
         results(k).u_non, ...
         '--', ...
         'LineWidth',line_width);

    yline(0,'k:','LineWidth',1.2);

    grid on;
    box on;

    xlabel('Time [s]', ...
        'FontSize',font_size);

    ylabel('Control Force [N]', ...
        'FontSize',font_size);

    title(sprintf('Initial angle = %d deg', ...
        initial_angles_deg(k)), ...
        'FontSize',title_size);

    legend('Linear','Nonlinear', ...
        'Location','best', ...
        'FontSize',10);

    set(gca, ...
        'FontSize',font_size, ...
        'LineWidth',1.0);

end

sgtitle('Linear vs Nonlinear Cart-Pole: LQR Control Effort', ...
    'FontSize',15, ...
    'FontWeight','bold');

%% ============================================================
% 12. Save figures as PDF
% ============================================================

exportgraphics(fig1, ...
    'Fig1_Pendulum_Angle.pdf', ...
    'ContentType','vector');

exportgraphics(fig2, ...
    'Fig2_Cart_Position.pdf', ...
    'ContentType','vector');

exportgraphics(fig3, ...
    'Fig3_Control_Effort.pdf', ...
    'ContentType','vector');

disp(' ');
disp('==============================================');
disp('Figures successfully saved as PDF');
disp('==============================================');

disp('Fig1_Pendulum_Angle.pdf');
disp('Fig2_Cart_Position.pdf');
disp('Fig3_Control_Effort.pdf');


%% ============================================================
% 13. Quantitative comparison
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('        NONLINEAR CART-POLE PERFORMANCE\n');
fprintf('============================================================\n');

fprintf('%10s %15s %15s %15s\n', ...
    'Angle', ...
    'Final Angle', ...
    'Max Angle', ...
    'Max Force');

fprintf('------------------------------------------------------------\n');

for k = 1:length(initial_angles_deg)

    theta = results(k).x_non(:,3);

    final_angle = abs(theta(end))*180/pi;

    max_angle = max(abs(theta))*180/pi;

    max_force = max(abs(results(k).u_non));

    fprintf('%10d %15.4f %15.4f %15.4f\n', ...
        initial_angles_deg(k), ...
        final_angle, ...
        max_angle, ...
        max_force);

end

fprintf('============================================================\n');


%% ============================================================
% LINEAR CART-POLE
% ============================================================

function dx = linear_cartpole(~,x,A,B,K)

    u = -K*x;

    dx = A*x + B*u;

end


%% ============================================================
% NONLINEAR CART-POLE
% ============================================================

function dx = nonlinear_cartpole(~,x,K,M,m,l,b,g)

    % States

    pdot  = x(2);
    theta = x(3);
    omega = x(4);

    %% LQR controller

    u = -K*x;

    %% Nonlinear dynamics

    s = sin(theta);
    c = cos(theta);

    denominator = M + m*s^2;

    %% Cart acceleration

    pddot = ...
        (u ...
        + m*s*(l*omega^2 + g*c) ...
        - b*pdot) ...
        / denominator;

    %% Angular acceleration

    theta_ddot = ...
        (-u*c ...
        - m*l*omega^2*c*s ...
        - (M+m)*g*s ...
        + b*pdot*c) ...
        / (l*denominator);

    %% State derivative

    dx = [pdot;
          pddot;
          omega;
          theta_ddot];

end