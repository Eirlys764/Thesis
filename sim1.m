close all;
clc;
clear all;

ftsz = 40;

s = @(angle) sin(angle);
c = @(angle) cos(angle);
Rz = @(angle) [c(angle) -s(angle); s(angle) c(angle)];

gamma = pi/3;
vM=5;

%% initial positions of relay and agent & Bearing Control gain
r0 = [0;0];
gbi = [0;1]-r0;
gbi = gbi/norm(gbi);
gFoV1 = Rz(gamma)*gbi; % negative x-axis
gFoV2 = Rz(-gamma)*gbi; % positive x-axis

x1_0 = [1,30]';
epsl = 2;

v1 = [0 -vM]';  % the agent moves along negative y-axis to escape

Kr_star = vM / s(gamma);
Kr = 1.5*Kr_star;
     

%% Gain selection for PI controller
xi = 0.5; % threshold of the steady-state error
alpha = 1; % must be >= 1   %Final error: error = 0.0004 
UB = ((1+xi)*epsl) / s(gamma); % pi_y_star(inf) - pr_y(inf)
Kp = alpha * vM*s(gamma)/(min(1,xi)*epsl); 
settling = 10;
Ki = 1*(Kp^2)/settling;    % don't choose settling = 4, why?

%% extra inertia component for brake action
lambda1 = -Kp/2 + sqrt((Kp^2) / 4 - Ki);
lambda2 = -Kp/2 - sqrt((Kp^2) / 4 - Ki);
Dr = max(epsl/s(gamma), norm(r0-x1_0));  % = e(0)
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

%% Solver
param.n =1;
param.v1 = v1;
param.gbi =gbi ;
param.Kr = Kr;
param.vM =vM ;
param.epsl =epsl ;
param.Kp = Kp;
param.Ki = Ki;
param.gFoV1 = gFoV1;
param.gFoV2 = gFoV2;
param.gamma = gamma;
param.BRAKE = BRAKE;

dt =10e-3;
tspan = (0:dt:30);

e_int_0 = zeros(2,1);
% e0 = err(r0, x1_0, epsl, gamma);
y0 = [x1_0; r0; e_int_0];

[t, y] = ode45(@(t,y) dyn1agent(t,y,param), tspan, y0);

x1 = y(:, 1:2); % arcoss the time
r = y(:, 3:4); 
e = compute_err_traj(r, x1, epsl, gamma);  % not from solver e=y(:, 7:8)
tt = t;

%% Count Violation (FoV constraint check)
% violation = 0;
% for kt = 1:length(t)
%     gr1 = x1(:, kt)-r(:, kt);
%     gr1 = gr1/norm(gr1);
%     check = (Rz(90)*gFoV1)'*gr1 <= 0 && (Rz(90)*gFoV2)'*gr1 >= 0;
%     if check == 0
%         fprintf("FoV constraint has been violated at t = %.4f seconds\n", t(kt));
%         violation = violation + 1;
%     end
% end


%% plot
% --- Plot 1: 2D Trajectories ---
figure
grid on
hold on
lw = 2.5;
plot(x1(:,1),x1(:,2),'b:','linewidth',lw/10)
plot(r(:,1),r(:,2),'r--','linewidth',lw/10)
set(gca,'FontSize',ftsz);

L = 5;
DeltaT = (tspan(end)-tspan(1))/L;
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
    0, 'color','blue','linewidth',lw,'MaxHeadSize',0.3);
    
    hx1 = plot(x1(kt,1),x1(kt,2), '-s', 'MarkerSize', 10, 'MarkerEdgeColor','blue', 'MarkerFaceColor','blue');
    hr  = plot(r(kt,1),r(kt,2),  '-s', 'MarkerSize', 10, 'MarkerEdgeColor','red',  'MarkerFaceColor','red');
    text(x1(kt,1)+0.05, x1(kt,2)+0.15, strcat('$k=',num2str(kl),'$'), ...
    'interpreter','latex','fontsize',ftsz*0.45,'Color','blue');
end


axis equal;
xlim([-10+min(x1(1,:)) 50]);
ylim([-30 10+max(r(2,:))]);
% xlim([min([x1(:,1); r(:,1)])-4, max([x1(:,1); r(:,1)])+4]);
% ylim([min([x1(:,2); r(:,2)])-1, max([x1(:,2); r(:,2)])+1]);
set(gca,'TickLabelInterpreter','latex', 'FontSize',ftsz);
hplots = [hr hx1 hgr1 hfov1];
leg = legend(hplots,{...
    '$\mathbf{p}_{r}(t_{k})$',...
    '$\mathbf{p}_{1}(t_{k})$',...
    '$\mathbf{g}_{r1}(t_{k})$',...
    '$\mathbf{h}_{FoVj}(t_{k}), ~ j=1,2$',...
    },'Interpreter','latex','location','northeast');
xlabel('$x$ [m]','interpreter','latex','fontsize',ftsz);
ylabel('$y$ [m]','interpreter','latex','fontsize',ftsz);
leg.FontSize = ftsz * 0.78;





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
function e = compute_err_traj(r, x1, epsl, gamma)

e = zeros(size(r));
for t =1:size(r,1)  
    e(t, :) = err(r(t, :)', x1(t, :)', epsl, gamma)';
end
end
