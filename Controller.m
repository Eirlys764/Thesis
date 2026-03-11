% p_r      : [2x1] relay position [x; y]
% p_i : [2 x n] matrix of agent positions
% gFOV1    : [2x1] direction vector of first FOV boundary
% gFOV2    : [2x1] direction vector of second FOV boundary

close all;
p_r0 = [0; 0]; %column vector
p_i = [1 2 3 4;
       2 4 6 8];   % agents are static
gamma = 30; % degree
gFOV1 = [cosd(90+gamma); sind(90+gamma)];
gFOV2 = [cosd(90-gamma); sind(90-gamma)];
epsilon = 0.05;
Kp = 5;

%% Find i^*, j^*
[i_star, j_star, ~] = compute_i_star(p_r0, p_i, gFOV1, gFOV2);

%% Find pr,y^*(t)
p_r_y_star = compute_p_r_y_star(p_r0, p_i, i_star, j_star, epsilon, gamma);

%% Move the relay to the desired position by minimizing the error to zero  
% Assume i^*=i^*(t), j^*=j^*(t)
% Assume u(t) = v(t) = 0, agents are static
% if e(t)>0, relay moves upwards; 
% if e(t)=0, relay arrives the desired position
% pick Kp, which decides the speed of convergence

dt = 0.01;
duration = 10;
t = 0: dt: duration;

p_r_y = zeros(size(t));
e = zeros(size(t));
p_r_y(1) = p_r0(2);

for n=1:length(t)-1
    e(n) = p_r_y_star - p_r_y(n);
    C_r_y = Kp * e(n);
    dp_r_y = C_r_y;
    p_r_y(n+1) = p_r_y(n) + dp_r_y*dt;
end

figure;
subplot(2,1,1);
plot(t, p_r_y, 'k--', 'LineWidth', 1.5); hold on; % if missing "hold on", the second plot will replace the first, then legend will mismatch(2 labels but only 1 plot)
plot(t, p_r_y_star*ones(size(t)), 'r--', 'LineWidth', 1.5);
title('Relay Position p_{r,y}(t)');
xlabel('Time (s)');
ylabel('p_{r,y}');
legend('p_r_y', 'p_r_y_star');
grid on;

subplot(2,1,2);
plot(t, e, 'k--', 'LineWidth', 1.5);
title('Error e(t)');
xlabel('Time (s)');
ylabel('e(t)');
grid on;

fprintf('The closest agent index: i_star = %d\n', i_star);
fprintf('The closest FoV index: j_star = %d\n', j_star);
fprintf('Desired position: p_r_y_star = %.4f \n', p_r_y_star);
fprintf('Final error: error = %.4f \n', e(end));

idx = find(abs(e)<1e-3, 1, "first");  %find(condition, number of return value, the order(if have many satisfied values, "first"=return the first satisfied value))
fprintf('Error first reaches zero at t = %.4f seconds\n', t(idx));
