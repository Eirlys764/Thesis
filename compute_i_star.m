function [i_star, j_star, distances] = compute_i_star(p_r, p_i, gFOV1, gFOV2)
% p_r      : [2x1] relay position [x; y]
% p_i : [2 x n] matrix of agent positions
% gFOV1    : [2x1] direction vector of first FOV boundary
% gFOV2    : [2x1] direction vector of second FOV boundary

n=size(p_i, 2); 

% Note: assume p_r(0) = p_r
% line 1
m1 = gFOV1(2) / gFOV1(1);
a1 = m1; 
b1 = 1;
c1 = -(p_r(2)-m1*p_r(1));

% line 2
m2 = gFOV2(2) / gFOV2(1);
a2 = m2; 
b2 = 1;
c2 = -(p_r(2)-m2*p_r(1));

A = [a1, a2];
B = [b1, b2];
C = [c1, c2];

distances = zeros(n, 2);

for i=1:n
    px = p_i(1,i);
    py = p_i(2,i);

    for j=1:2
        distances(i, j)= abs(A(j)*px + B(j)*py + C(j))/sqrt(A(j)^2 + B(j)^2);
    end
end

% [the min value, the index] = min(A, B) A, B are comparable arrays, if B
% is empty(write as a placeholder []), add the third argument "dimension"
% that make A compares its own elements along that dimension, where '1' is
% column-wise, "2" is row-wise
[min_dist_per_agent, j_star_per_agent] = min(distances, [], 2);  % each row min[di1, di2]=di1, j= 1,2   
[~, i_star] = min(min_dist_per_agent); % min_dist_per_agent and j_star_per_agent are nx1 column vectors
j_star = j_star_per_agent(i_star); 