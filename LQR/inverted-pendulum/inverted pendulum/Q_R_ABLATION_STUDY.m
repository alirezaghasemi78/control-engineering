%% Q/R ABLATION STUDY
% Study the effect of Q and R on LQR performance

clear;
clc;
close all;

%% Plant

m = 1;
l = 1;
b = 0.2;
g = 9.81;

A = [0        1;
     g/l   -b/(m*l^2)];

B = [0;
     1/(m*l^2)];

%% Simulation settings

x0 = [10*pi/180; 0];

t = 0:0.001:5;

%% Baseline

Q_base = diag([100 10]);
R_base = 1;

[K_base,~,P_base] = lqr(A,B,Q_base,R_base);

%% ============================================================
% 1. Q sensitivity
% =============================================================

Q_scales = [0.1 0.5 1 2 10];

figure;

for i = 1:length(Q_scales)

    Q = Q_scales(i)*Q_base;

    K = lqr(A,B,Q,R_base);

    Acl = A-B*K;

    sys_cl = ss(Acl,B,eye(2),zeros(2,1));

    [x,t_sim] = initial(sys_cl,x0,t);

    plot(t_sim,x(:,1),'LineWidth',1.3);
    hold on;

end

grid on;

xlabel('Time [s]');
ylabel('\theta [rad]');

legend( ...
    '0.1Q',...
    '0.5Q',...
    '1Q',...
    '2Q',...
    '10Q');

title('Sensitivity to Q');

%% ============================================================
% 2. R sensitivity
% =============================================================

R_scales = [0.1 0.5 1 2 10];

figure;

for i = 1:length(R_scales)

    R = R_scales(i)*R_base;

    K = lqr(A,B,Q_base,R);

    Acl = A-B*K;

    sys_cl = ss(Acl,B,eye(2),zeros(2,1));

    [x,t_sim] = initial(sys_cl,x0,t);

    plot(t_sim,x(:,1),'LineWidth',1.3);
    hold on;

end

grid on;

xlabel('Time [s]');
ylabel('\theta [rad]');

legend( ...
    '0.1R',...
    '0.5R',...
    '1R',...
    '2R',...
    '10R');

title('Sensitivity to R');

%% ============================================================
% 3. Q effect on control effort
% =============================================================

figure;

for i = 1:length(Q_scales)

    Q = Q_scales(i)*Q_base;

    K = lqr(A,B,Q,R_base);

    Acl = A-B*K;

    sys_cl = ss(Acl,B,eye(2),zeros(2,1));

    [x,t_sim] = initial(sys_cl,x0,t);

    u = -(K*x')';

    plot(t_sim,u,'LineWidth',1.3);
    hold on;

end

grid on;

xlabel('Time [s]');
ylabel('u');

legend( ...
    '0.1Q',...
    '0.5Q',...
    '1Q',...
    '2Q',...
    '10Q');

title('Effect of Q on Control Effort');

%% ============================================================
% 4. R effect on control effort
% =============================================================

figure;

for i = 1:length(R_scales)

    R = R_scales(i)*R_base;

    K = lqr(A,B,Q_base,R);

    Acl = A-B*K;

    sys_cl = ss(Acl,B,eye(2),zeros(2,1));

    [x,t_sim] = initial(sys_cl,x0,t);

    u = -(K*x')';

    plot(t_sim,u,'LineWidth',1.3);
    hold on;

end

grid on;

xlabel('Time [s]');
ylabel('u');

legend( ...
    '0.1R',...
    '0.5R',...
    '1R',...
    '2R',...
    '10R');

title('Effect of R on Control Effort');


function metrics = evaluate_lqr(A,B,Q,R,x0,t)

    [K,~,P] = lqr(A,B,Q,R);

    Acl = A-B*K;

    sys_cl = ss(Acl,B,eye(size(A)),zeros(size(A,1),size(B,2)));

    [x,t] = initial(sys_cl,x0,t);

    u = -(K*x')';

    %% State cost

    state_cost_signal = sum((x*Q).*x,2);

    %% Control cost

    control_cost_signal = sum((u*R).*u,2);

    %% Total cost

    Jx = trapz(t,state_cost_signal);
    Ju = trapz(t,control_cost_signal);

    J = Jx + Ju;

    %% Metrics

    metrics.K = K;
    metrics.poles = P;

    metrics.Jx = Jx;
    metrics.Ju = Ju;
    metrics.J = J;

    metrics.max_state = max(abs(x(:,1)));

    metrics.rms_control = sqrt(mean(u.^2));

    metrics.max_control = max(abs(u));

    %% Settling time

    theta = abs(x(:,1));

    threshold = 0.02*abs(x0(1));

    idx = find(theta > threshold);

    if isempty(idx)

        metrics.settling_time = 0;

    else

        metrics.settling_time = t(idx(end));

    end

end

Q_scales = [0.1 0.5 1 2 10];

results = zeros(length(Q_scales),5);

for i = 1:length(Q_scales)

    Q = Q_scales(i)*Q_base;

    M = evaluate_lqr(A,B,Q,R_base,x0,t);

    results(i,:) = [ ...
        Q_scales(i), ...
        M.J, ...
        M.Jx, ...
        M.Ju, ...
        M.max_control];

end

disp("   Q scale    |    J |     Jx |    Ju |     max(u)")
disp(results);