function dydt = dyn1agent(t, y, param)

x1 = y(1: 2);  % the agent
r = y(3: 4); % the relay 
e_int = y(5:6); % the intergral 

%% tracjectory of the relay under PI controller and bearing-based controller

% Relay dynamics
gr1_ = x1-r;
gr1 = gr1_/norm(gr1_);

Pgr1 = eye(2) - gr1*gr1';
u_r = -param.Kr*Pgr1*param.gbi; 
[e, i_star, j_star] = err(r, x1, param.epsl, param.gamma); % by solving the diff of e, we get the intergral of e

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

% Single agent's dynamics
dotpi = [0 -param.vM]'; % assume it moves along y negative axis at vM
% dotpi = [0 0]'; % static

R = 1;
w = param.vM;
% dotpi = [-R*w*sin(w*t) 0]'; % harmonic motion: horizontal oscillation, y component=0
% dotpi = [-R*w*sin(w*t) R*w*cos(w*t)]'; % unite circular

dydt = [reshape(dotpi, 2*param.n, 1); dotpr; e_int_dot];  %reshape(A, #rows, #cols)
end

%% add at the end of the file, so everything after it wouldn't be treated as part of a new nested function 
function value = neg_pass(value)
if value > 0
    value = 0;
end

end



