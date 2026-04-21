function dydt = dynNagent(t, y, param)

n = param.n;
x = zeros(2,n);
gr = zeros(2,n);
r = y(2*n+1: 2*n+2);
e_int =y(2*n+3: 2*n+4);
for i=0: n-1
    ii = i+1;
    x(:,ii) = y(2*i+1: 2*i+2);
    gr(:,ii) = x(:,ii) - r;
    gr(:,ii) = gr(:,ii) / norm(gr(:,ii));
end

%% find the closest two agents to the FoV1,2
grigFov1 = -Inf;
grigFov2 = -Inf;
k1 =0;
k2 =0;

for i=1:n
    grigFov1_i = gr(:,i)'*param.gFoV1;
    grigFov2_i = gr(:,i)'*param.gFoV2;
    if grigFov1_i > grigFov1
        grigFov1 = grigFov1_i;
        k1 =i;
    end
    if grigFov2_i > grigFov2
        grigFov2 = grigFov2_i;
        k2 =i;
    end
end

k=k1;
if grigFov2 > grigFov1
   k=k2;
end


%% compute u_r = bearing-based controller
Xn = side(gr, param.gFoV1, param.gFoV2, n);
if Xn >= 0
    grbar = gr(:, k);
    Pgrbar = eye(2) - grbar*grbar';
    u_r = -param.Kr * Pgrbar * param.gbi;
else
    grbar_1 = gr(:, k1);
    grbar_2 = gr(:, k2);
    Pgrbar_1 = eye(2) - grbar_1*grbar_1';
    Pgrbar_2 = eye(2) - grbar_2*grbar_2';
    u_r = -param.Kr * (Pgrbar_1 + Pgrbar_2) * param.gbi;
end


%% PI controller
[e, i_star, j_star] = err(r, x, param.epsl, param.gamma); % by solving the diff of e, we get the intergral of e

Cr = param.Kp*e + param.Ki*e_int + param.BRAKE*[0, neg_pass(e(2))]';
dotpr_unsat = u_r + Cr;

% Anti-windup
thr = 3 * param.vM;
if norm(dotpr_unsat) > thr
    dotpr = thr * dotpr_unsat / norm(dotpr_unsat);
    e_int_dot = zeros(2,1);  % freeze integrator
else
    dotpr = dotpr_unsat;
    e_int_dot = e;            % normal integration
end


%% n agents dynamic
R = param.R;
w = param.vM;
% dotpi_1 = [0 -param.vM]'; % moves along negative y-axis
dotpi_1 = [0 0]';
dotpi_2 = [-R*w*sin(w*t) 0]'; % harmonic motion: horizontal oscillation, y component=0
dotpi_3 = [-R*w*sin(w*t)*exp(-2*t) 0]'; 
% dotpi_2 = [t -R*w*sin(w*t)]'; 
% dotpi_3 = [t -R*w*sin(w*t)*exp(-2*t)]'; 
dotpi_4 = [-R*w*sin(w*t) R*w*cos(w*t)]'; % unite circular
dotpi_5 = [-R*w*sin(w*t)*exp(-2*t) R*w*cos(w*t)*exp(-2*t)]'; 
% dotpi = [dotpi_1, dotpi_2, dotpi_3, dotpi_4, dotpi_5];


%% differential components
dydt = [dotpi_1; dotpi_2; dotpi_3; dotpi_4; dotpi_5; dotpr; e_int_dot];  
end
% fprintf('The closest agent index: i_star = %d\n', i_star);
% fprintf('The closest FoV index: j_star = %d\n', j_star);

function value = neg_pass(value)
if value > 0
    value = 0;
end
end

