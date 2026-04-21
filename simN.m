% close all;
% clc;
clear all;

ftsz = 40;

s = @(angle) sin(angle);
c = @(angle) cos(angle);
Rz = @(angle) [c(angle) -s(angle); s(angle) c(angle)];

gamma = pi/3;
vM=5;
epsl = 2;
R=1;

% initial positions of relay and agent
r0 = [0;-10];
gbi = [0;1]-r0;
gbi = gbi/norm(gbi);
gFoV1 = Rz(gamma)*gbi; % negative x-axis
gFoV2 = Rz(-gamma)*gbi; % positive x-axis
x1_0 = [R 3]'; % pi_y_star = pr(0)
x2_0 = [2*R,3]';
x3_0 = x2_0;
x4_0 = [-2*R,3]';
x5_0 = x4_0;

% x1_0 = r0+20*gFoV1; % agent 1 is initially located on first border
% x2_0 = r0+27*Rz(pi/100)*gFoV1;
% x3_0 = r0+20*Rz(pi/8)*gFoV1; % chattering dance
% x4_0 = r0+6*Rz(-pi/6)*gbi; % circle
% x5_0 = r0+30*Rz(pi/8)*gFoV1;

% temp_pos = r0+30*Rz(pi/8)*gFoV1;
% x5_0 = [temp_pos(1); -10];  % Lorenz strange attractor

x_0 = [x1_0, x2_0, x3_0, x4_0, x5_0];


%%  Gain selection for bearing-based controller
g_gamma_star = 0;
if gamma >0 && gamma <= pi/6
    g_gamma_star = 2*s(gamma)^3;
elseif gamma >pi/6 && gamma <= pi/2
    g_gamma_star = (3*s(gamma)-1) / 2;
end

Kr_star = vM / g_gamma_star;
Kr = 1.5*Kr_star;

%% Gain selection for PI controller
xi = 0.5; % threshold of the steady-state error
alpha = 20; % must be >= 1  %Final error: error = 0.0197
UB = ((1+xi)*epsl) / s(gamma); % pi_y_star(inf) - pr_y(inf)
Kp = alpha * vM*s(gamma)/(min(1,xi)*epsl); 
settling = 10;
Ki = 1*(Kp^2)/settling;


%% extra inertia component for brake action
lambda1 = -Kp/2 + sqrt((Kp^2) / 4 - Ki);
lambda2 = -Kp/2 - sqrt((Kp^2) / 4 - Ki);

Dr =0;
for i=1:size(x_0, 2) % #cols of x_0
    Dr = max(epsl/s(gamma), norm(r0-x_0(:,i))); % x_0(:,i) extra each col=vector (x1,y1)
end

z = (vM-lambda2*Dr) / (vM-lambda1*Dr);
z1 = z^(lambda1 / (lambda1 - lambda2));
z2 = z^(lambda2 / (lambda1 - lambda2));
Emax = (z1*(Dr-vM/lambda1) - z2*(Dr-vM/lambda2)) / (lambda1-lambda2) - vM/(lambda1*lambda2);

if abs(imag(Emax)) < 1e-6
    Emax = real(Emax);
end

BRAKE = Ki*Emax*s(gamma)/epsl;
if BRAKE < 0
    error("BRAKE is wrong");
end
if Ki == 0
    BRAKE =0;
end


%% ode45
param.n =5;
param.gbi =gbi ;
param.Kr = Kr;
param.vM =vM ;
param.epsl =epsl ;
param.Kp = Kp;
param.Ki = Ki;
param.gFoV1 = gFoV1;
param.gFoV2 = gFoV2;
param.gamma = gamma;
param.R = R;
param.BRAKE = BRAKE;

dt =10e-3;
tspan = (0:dt:30);

e_int_0 = zeros(2,1);
y0 = [x1_0; x2_0; x3_0; x4_0; x5_0; r0; e_int_0];
[t, y] = ode45(@(t,y) dynNagent(t,y,param), tspan, y0);

x1 = y(:, 1:2); % arcoss the time % 3001x2 
x2 = y(:, 3:4);
x3 = y(:,5:6);
x4 = y(:,7:8);
x5 = y(:,9:10);
r = y(:,11:12);
e = compute_err_traj(r, x1,x2,x3,x4,x5, epsl, gamma);
tt = t; 


%% plot
% --- Plot 1: 2D Trajectories ---
figure
grid on
hold on
lw = 2.5;
plot(x1(:,1),x1(:,2),'b:','linewidth',lw/10)
plot(x2(:,1),x2(:,2),'m:','linewidth',lw/10)
plot(x3(:,1),x3(:,2),'color','yellow','linewidth',lw/10)
plot(x4(:,1),x4(:,2),'color','cyan','linewidth',lw/10)
plot(x5(:,1),x5(:,2),'color','green','linewidth',lw/10)
plot(r(:,1),r(:,2),'r--','linewidth',lw/10)
set(gca,'FontSize',ftsz);

L = 5;
DeltaT = (tt(end)-tt(1))/L;
hc = 150;
ls = linspace(0,1);
seg1 = zeros(2,length(ls));
seg2 = zeros(2,length(ls));
for kl = 0:L
    kt = round(kl*DeltaT/dt);
    if kt < 1
        kt = 1;
    end
  
    x1_kt = x1(kt,:)';  % 2×1
    x2_kt = x2(kt,:)';
    x3_kt = x3(kt,:)';
    x4_kt = x4(kt,:)';
    x5_kt = x5(kt,:)';
    r_kt  = r(kt,:)';   % 2×1 
    fov1_kt = r_kt+gFoV1*hc;
    fov2_kt = r_kt+gFoV2*hc;


    for kk = 1:length(ls)
        seg1(:,kk) = ls(kk)*r_kt+(1-ls(kk))*fov1_kt;
        seg2(:,kk) = ls(kk)*r_kt+(1-ls(kk))*fov2_kt;
    end
    hfov1 = plot(seg1(1,:),seg1(2,:),':k','linewidth',lw/1.5);
    hfov2 = plot(seg2(1,:),seg2(2,:),':k','linewidth',lw/1.5);
    hgr1 = quiver(r_kt(1), r_kt(2),x1_kt(1)-r_kt(1), x1_kt(2)-r_kt(2), ...
    0, 'color','blue','linewidth',lw,'MaxHeadSize',0.2);
    hgr2 = quiver(r_kt(1), r_kt(2),x2_kt(1)-r_kt(1), x2_kt(2)-r_kt(2), ...
    0, 'color','magenta','linewidth',lw,'MaxHeadSize',0.2);
    hgr3 = quiver(r_kt(1), r_kt(2),x3_kt(1)-r_kt(1), x3_kt(2)-r_kt(2), ...
    0, 'color','yellow','linewidth',lw,'MaxHeadSize',0.2);
    hgr4 = quiver(r_kt(1), r_kt(2),x4_kt(1)-r_kt(1), x4_kt(2)-r_kt(2), ...
    0, 'color','cyan','linewidth',lw,'MaxHeadSize',0.2);
    hgr5 = quiver(r_kt(1), r_kt(2),x5_kt(1)-r_kt(1), x5_kt(2)-r_kt(2), ...
    0, 'color','green','linewidth',lw,'MaxHeadSize',0.2);
    
    hx1 = plot(x1(kt,1),x1(kt,2), '-s', 'MarkerSize', 10, 'MarkerEdgeColor','blue', 'MarkerFaceColor','blue');
    hx2 = plot(x2(kt,1),x2(kt,2), '-s', 'MarkerSize', 10, 'MarkerEdgeColor','magenta', 'MarkerFaceColor','magenta');
    hx3 = plot(x3(kt,1),x3(kt,2), '-s', 'MarkerSize', 10, 'MarkerEdgeColor','yellow', 'MarkerFaceColor','yellow');
    hx4 = plot(x4(kt,1),x4(kt,2), '-s', 'MarkerSize', 10, 'MarkerEdgeColor','cyan', 'MarkerFaceColor','cyan');
    hx5 = plot(x5(kt,1),x5(kt,2), '-s', 'MarkerSize', 10, 'MarkerEdgeColor','green', 'MarkerFaceColor','green');
    hr  = plot(r(kt,1),r(kt,2),  '-s', 'MarkerSize', 10, 'MarkerEdgeColor','red',  'MarkerFaceColor','red');
    text(r(kt,1)+0.05, r(kt,2)+0.15, strcat('$k=',num2str(kl),'$'), ...
    'interpreter','latex','fontsize',ftsz*0.45,'Color','blue');
end


axis equal;
% xlim([-10+min(x1(1,:)) 50]);
% ylim([-30 10+max(r(2,:))]);
x = [x1; x2; x3; x4; x5];
xlim([min([x(:,1);  r(:,1)])-3, max([x(:,1); r(:,1)])+3]);
ylim([min([x(:,2);  r(:,2)])-1, max([x(:,2); r(:,2)])+1]);
set(gca,'TickLabelInterpreter','latex', 'FontSize',ftsz);
hplots = [hr hx1 hx2 hx3 hx4 hx5 hgr1 hgr2 hgr3 hgr4 hgr5 hfov1];
leg = legend(hplots,{...
    '$\mathbf{p}_{r}(t_{k})$',...
    '$\mathbf{p}_{1}(t_{k})$',...
    '$\mathbf{p}_{2}(t_{k})$',...
    '$\mathbf{p}_{3}(t_{k})$',...
    '$\mathbf{p}_{4}(t_{k})$',...
    '$\mathbf{p}_{5}(t_{k})$',...
    '$\mathbf{g}_{r1}(t_{k})$',...
    '$\mathbf{g}_{r2}(t_{k})$',...
    '$\mathbf{g}_{r3}(t_{k})$',...
    '$\mathbf{g}_{r4}(t_{k})$',...
    '$\mathbf{g}_{r5}(t_{k})$',...
    '$\mathbf{h}_{FoVj}(t_{k}), ~j=1,2$',...
    },'Interpreter','latex','location','northeast');
xlabel('$x$ [m]','interpreter','latex','fontsize',ftsz);
ylabel('$y$ [m]','interpreter','latex','fontsize',ftsz);
leg.FontSize = ftsz * 0.4;





% --- Plot 2: Error over time ---
figure;
hplots = plot(tt,e);
set(gca,'TickLabelInterpreter','latex', 'FontSize',ftsz);
leg = legend(hplots,{...
    '$e_x(t)$',...
    '$e_y(t)$',...
    },'Interpreter','latex','location','southeast');
xlabel('$t$ [s]','interpreter','latex','fontsize',ftsz);
ylabel('$e(t)$ [m]','interpreter','latex','fontsize',ftsz);


fprintf('Final error: error = %.4f \n', e(end));

%%
function e = compute_err_traj(r, x1, x2,x3,x4,x5, epsl, gamma)

e = zeros(size(r));
for t =1:size(r,1)
    e(t, :) = err(r(t, :)', [x1(t, :)', x2(t, :)',x3(t, :)',x4(t, :)',x5(t, :)'], epsl, gamma)';
end
end