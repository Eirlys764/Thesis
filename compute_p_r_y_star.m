function p_r_y_star = compute_p_r_y_star(p_r, p_i, i_star, j_star, epsilon, gamma)
% epsilon  : min safety distance 
% gamma    : half-angle of the FOV (in radians)
% p_r      : [2x1] relay position [x; y]
% p_i : [2 x n] matrix of agent positions

% p_r_x = p_r(1);
p_r_y = p_r(2);

p_i_star_x = p_i(1, i_star);
p_i_star_y = p_i(2, i_star);

cot(gamma) = cos(gamma)/sin(gamma);

p_r_y_star = p_r_y - epsilon/sin(gamma) + ((-1)^j_star)*cot(gamma)*p_i_star_x + p_i_star_y;
end