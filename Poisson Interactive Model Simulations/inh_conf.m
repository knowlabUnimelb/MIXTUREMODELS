% NOTE: high salient distractors are leading to faster completion times
% than low salient distractors - why?

function cdf = inh_conf(params, varargin)
% Poisson parallel interactive model with inhibitory connections and an
% assumption of conflict
%
% The poisson distribution expresses the probability of a given number of 
%  events occurring in a fixed interval of time and/or space if these events 
%  occur with a known average rate and independently of the time since the 
%  last event
%
% This function implements the transition matrix method outlined in Table A.1 
% of Johnson, Blaha, Houpt & Townsend (2010) JMP, 54, 53-72

%% Optional arguments
optargs = {(0:1000)};
newVals = cellfun(@(x) ~isempty(x), varargin);% skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
[t] = optargs{:}; % Place optional args in memorable variable names

%%
% h  :  Time scaling constant (should be small)
v_a  = params(1); % v_a:  Rate in channel A ( multiplied by time scale )
v_d  = params(2); % v_b:  Rate in channel B ( multiplied by time scale )
p_ad = params(3); % p  :  Probability of sharing from a to b
p_da = params(4); % p  :  Probability of sharing from b to a
c_a  = params(5); % A threshold
c_d  = params(6); % B threshold

%% Set up transition matrix
% The logic is that a channel can generate an increment with some
% probability, that increment also has some probability of being shared. If
% shared, the other channel decrements. If a channel generates an
% increment, you need to specify what happens to it (shared or not shared).
% Sharing can only occur when an increment is generated. (i.e., can't go
% from 0 to -1)

% P[0,0] = P(A & B Same)
p00 = (1-v_a) * (1-v_d) +...
      v_a  * v_d * p_ad * p_da;
p00_ac = 1-v_d; % if a has reached criterion, then response will terminate
p00_dc = 1-v_a; % if d has reached criterion, then nothing will happen

% P[1,0] = P(A Up & B Same)
p10 = v_a * v_d * p_ad * (1-p_da) +... % A increments, D increments, A inhibits D, D doesn't inhibit A
      v_a * (1-v_d) * (1-p_ad);        % or A increments, D doesn't increment, A doesn't inhibit D
p10_d0 = v_a * (1-v_d) * p_ad; % p1_1 % added if d has no counts: A increments, D doesn't increment, A inhibits D (same as p1_1)
p10_dc = v_a; % if d has reached criterion

% P[0,1] = P(A Same & B Up)
p01 = v_a * v_d * p_da * (1-p_ad) +... % A increments, D increments, D inhibits A, A doesn't inhibit D
     (1-v_a) * v_d * (1-p_da);         % or A doesn't increment, D increments, and D doesn't inhibit A
p01_a0 = (1-v_a) * v_d * p_da; % p_11 % added if a has no counts: A doesn't increment, D increments, D inhibits A
p01_ac = v_d; % if a has reached criterion 

% P[1,1] = P(A Up & B Up)
p11 = v_a * v_d * (1-p_ad) * (1-p_da); % A and D increment and nobody inhibits the other

% P[0,-1] = P(A Same & B Down)
p0_1 = 0;

% P[-1,0] = P(A Down & B Same)
p_10 = 0;

% P[1,-1] = P(A Up & B Down)
p1_1 = v_a * (1-v_d) * p_ad; % A increments, D doesn't increment, A inhibits D

% P[-1,1] = P( A Down & B Up)
p_11 = (1-v_a) * v_d * p_da; % A doesn't increment, D increments, D inhibits A

% P[-1,-1] = P(A Same & B Down)
p_1_1 = 0;


%%
% Main diagonal component
T1_0  = diag([p00 * ones(c_d,1); p00 + p01]) +... % Diagonal component
         diag(p01 * ones(c_d,1), 1);              % right shift one diagonal component
% Main diagonal component if a state = 0     
T1_a0 = diag([p00 * ones(c_d,1); p00+p01+p_11]) +... % Diagonal component
        diag((p01+p_11) * ones(c_d,1), 1);            % right shift one diagonal component
% Main diagonal component if a state = c_a    
T1_ac = diag([ones(c_d+1,1)]);                         % 1 matrix

% Right shift one component
T2_0  = diag([p1_1*ones(c_d,1)], -1) + ...                  % Left shift one diagonal
        diag([p10+p1_1; p10*ones(c_d-1,1); p10+ p11]) + ... % Main diagonal
        diag(p11 * ones(c_d,1), 1);                         % Right shift one diagonal
     
% Left shift one component    
T3_0  = diag([zeros(c_d,1); p_11]) + ...  % Main diagonal
        diag(p_11 * ones(c_d,1), 1);      % Left shift one diagonal

T1 = mat2cell(repmat(T1_0, 1, c_a-1), c_d+1, repmat(c_d+1, 1, c_a-1));
T2 = mat2cell(repmat(T2_0, 1, c_a-1), c_d+1, repmat(c_d+1, 1, c_a-1));
% T3 = mat2cell(repmat(T3_0, 1, c_a-2), c_b+1, repmat(c_b+1, 1, c_a-2));
T3 = mat2cell(repmat(T3_0, 1, max(1,c_a-2)), c_d+1, repmat(c_d+1, 1, max(1,c_a-2)));

%% Build full transition matrix
T = blkdiag(T1_a0, T1{1:c_a-1}, T1_ac) + ...                        
    vertcat(blkdiag(horzcat(zeros(c_d+1,c_d+1), T2_0), T2{1:c_a-1}), ...
       zeros(c_d+1,(c_a+1)*(c_d+1))) + ...
    horzcat(blkdiag(vertcat(zeros(c_d+1,c_d+1), T3{1}), T3{1:c_a-2}, ...
       zeros(c_d+1,c_d+1)), zeros((c_d+1)*(c_a+1),c_d+1));
   
%% Internal transitions
k = size(T,1);             % Size of matrix
% Q = T(1:(k-1), 1:(k-1));   % Transient states
Q = T(1:(c_a+1)*c_d, 1:(c_a+1)*c_d);

%% Absorbing states
R = T(1:(c_a+1)*c_d, ((c_a+1)*c_d+1):k);         % Absorbing state vector

%%
k = size(Q,1);             % Size of transition matrix
cdf = zeros(size(t));      % Vector of probilities at each time point
I = eye(k);

%% Compute cdf
% The more general form allows for a starting point vector of Z. Here we
% assume that accumulation starts from the first position (0, 0)

Q_mult = Q^(t(1) + 1);                      % Q^(N+1), N = t(1)
C_N = sum(inv(I - Q) * (I-Q_mult) * R, 2);  % (I - Q)^-1 * (I - Q^(N+1)) * R
cdf(1) = C_N(1);                           % CDF

for i = 2:size(t,2)
    Q_mult = Q_mult * Q^(t(i) - t(i-1));         % Q^(N+1)
    C_N = sum(inv(I - Q) * (I - Q_mult) * R, 2);
    cdf(i) = C_N(1);
end