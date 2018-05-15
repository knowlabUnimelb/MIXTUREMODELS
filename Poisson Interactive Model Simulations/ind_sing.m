function C_Ni = ind_sing(params, t_vec)

v_a  = params(1);
c_a  = params(2);


C_Ni = 1-poisscdf(c_a-1, v_a*t_vec); % CDF
% C_Ni = poisscdf(c_a-1, v_a*t_vec); % Survivor