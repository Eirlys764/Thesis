function dydt = dyn2agent(t, y, param)

x1 = y(1:2);
x2 = y(3:4);
x = [x1 x2];
r = y(5:6);
e_int = y(7:8);

gr = [(x1-r)/norm(x1-r) (x2-r)/norm(x2-r)]; %2x2 matrix
gr1 = gr(:, 1);
gr2 = gr(:, 2);

Pgbi = eye(2)-param.gbi*param.gbi';
Pgr1 = eye(2)-gr1*gr1';
Pgr2 = eye(2)-gr2*gr2';

X = sign(gr1' * Pgbi*gr2);

%% bearing-based controller
if X>=0
    gr1gFovj = max(gr1'*param.gFoV1, gr1'*param.gFoV2);
    gr2gFovj = max(gr2'*param.gFoV1, gr2'*param.gFoV2); % the two should on the same side!
    grbar = gr2;
    if gr1gFovj > gr2gFovj
        grbar = gr1;
    end
    
    Pgrbar = eye(2)-grbar*grbar';
    u_r = -param.Kr* Pgrbar * param.gbi;

else
    u_r = -param.Kr* (Pgr1 + Pgr2) * param.gbi;
end

%% PI controller
[e, i_star, j_star] = err(r, x, param.epsl, param.gamma); % by solving the diff of e, we get the intergral of e

%% 2 agents dynamic
R = param.R;
w = param.vM;
dotpi_1 = [-R*w*sin(w*t)*exp(-2*t) 0]'; % harmonic motion: horizontal oscillation, y component=0
dotpi_2 = [-R*w*sin(w*t) R*w*cos(w*t)]'; % unite circular

%% differential components
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

dydt = [dotpi_1; dotpi_2; dotpr; e_int_dot];  %reshape(A, #rows, #cols)
end

%%
function value = neg_pass(value)
if value > 0
    value = 0;
end

end