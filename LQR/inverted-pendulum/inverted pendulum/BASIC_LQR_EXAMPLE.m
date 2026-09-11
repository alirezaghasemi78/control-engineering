%% BASIC LQR EXAMPLE
% LQR control of an unstable inverted pendulum
% ------------------------------------------------

clear;
clc;
close all;

%% 1. Plant parameters

m = 1;          % Mass [kg]
l = 1;          % Pendulum length [m]
b = 0.2;        % Damping coefficient
g = 9.81;       % Gravity [m/s^2]

%% 2. State-space model

A = [0        1;
     g/l   -b/(m*l^2)];

B = [0;
     1/(m*l^2)];

C = eye(2);
D = zeros(2,1);

%% 3. Check controllability

Co = ctrb(A,B);

fprintf('Rank of controllability matrix = %d\n',rank(Co));
fprintf('Number of states = %d\n',size(A,1));

%% 4. Open-loop poles

p_open = eig(A);

fprintf('\nOpen-loop poles:\n');
disp(p_open);

%% 5. LQR weights

Q = diag([100 10]);
R = 1;

%% 6. LQR design

[K,S,P] = lqr(A,B,Q,R);

fprintf('\nLQR gain K:\n');
disp(K);

fprintf('Riccati solution S:\n');
disp(S);

fprintf('Closed-loop poles:\n');
disp(P);

%% 7. Closed-loop system

Acl = A - B*K;

sys_open = ss(A,B,C,D);
sys_closed = ss(Acl,B,C,D);

%% 8. Initial condition

x0 = [10*pi/180; 0];

%% 9. Time vector

t = 0:0.001:5;

%% 10. Open-loop response

[y_open,t_open,x_open] = initial(sys_open,x0,t);

%% 11. Closed-loop response

[y_closed,t_closed,x_closed] = initial(sys_closed,x0,t);

%% 12. Control input

u_closed = -(K*x_closed')';


%% 13. پاسخ زاویه (شکل شماره 1) - تفکیک شده در دو زیرنمودار
figure;

subplot(2,1,1);
plot(t_open, x_open(:,1), 'b', 'LineWidth', 1.5);
grid on;
title('Open-loop Response (Unstable)');
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');

subplot(2,1,2);
plot(t_closed, x_closed(:,1), 'b', 'LineWidth', 1.5);
grid on;
title('Closed-loop Response (LQR)');
xlabel('Time [s]');
ylabel('$\theta$ [rad]', 'Interpreter', 'latex');

sgtitle('Figure 1: Angle Response (\theta)');

%% 14. پاسخ سرعت زاویه‌ای (شکل شماره 2) - تفکیک شده در دو زیرنمودار
figure;

subplot(2,1,1);
plot(t_open, x_open(:,2), 'r', 'LineWidth', 1.5);
grid on;
title('Open-loop Response (Unstable)');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

subplot(2,1,2);
plot(t_closed, x_closed(:,2), 'r', 'LineWidth', 1.5);
grid on;
title('Closed-loop Response (LQR)');
xlabel('Time [s]');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

sgtitle('Figure 2: Angular Velocity Response (\dot{\theta})');

%% 15. Control effort (شکل شماره 3) - فقط اصلاح خطای Interpreter
figure;

plot(t_closed, u_closed, 'k', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Control input $u$', 'Interpreter', 'latex');  % رفع خطا
title('Figure 3: LQR Control Effort');

%% 16. Phase portrait (شکل شماره 4) - تفکیک شده در دو زیرنمودار
figure;

subplot(1,2,1);
plot(x_open(:,1), x_open(:,2), 'b', 'LineWidth', 1.5);
grid on;
axis equal;  % نسبت ابعاد یکسان برای نمایش صحیح دایره‌وار
title('Open-loop Phase');
xlabel('$\theta$ [rad]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

subplot(1,2,2);
plot(x_closed(:,1), x_closed(:,2), 'r', 'LineWidth', 1.5);
grid on;
axis equal;
title('LQR Closed-loop Phase');
xlabel('$\theta$ [rad]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex');

sgtitle('Figure 4: Phase Portrait Comparison');

%% 17. Cost calculation

state_cost = sum(sum((x_closed * Q) .* x_closed,2));
control_cost = sum(u_closed.^2 * R);

J = trapz(t_closed, ...
    sum((x_closed * Q).*x_closed,2) + ...
    (u_closed.^2)*R);

fprintf('\nApproximate LQR cost J = %.4f\n',J);








