function [e, i_star, j_star] = err(p_r, p_i, epsilon, gamma)
% This function computes the error between the reference signal and the 
% actual position of the relay drone.

% p_i: 2xn nagents matrix

w_plus = [1/tan(gamma), -1];
w_minus = [-1/tan(gamma), -1];

j_star =0;
i_star =0;
term = +Inf;

% at the end of the optimization we expect that the value of term is: 
% term = <[(-1)^j_star / tan(gamma), -1], (p_r-p_i_star)>
for i=1:size(p_i,2)
    p_r_p_i = p_r - p_i(:,i); % select all the rows of the ith column
    [term_ij, j] = min([abs(w_minus*p_r_p_i) abs(w_plus*p_r_p_i)]); %j=1 or 2 
    % [a b] is the syntax for horizontal concatenation.
    % min([a b]): passing one input argument: a vector containing two values
    % If pass two separate arguments, MATLAB performs an element-wise comparison.
    % If A and B are matrices, it returns a new matrix of the same size where each element is the smaller of the two corresponding elements.
    if term_ij < term
        term = term_ij;
        j_star = j;
        i_star = i;
    end
end

e = [0, -epsilon/sin(gamma)+term]';
end