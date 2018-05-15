function cdf = inh_and(params, varargin)
% Poisson parallel interactive model with inhibitory connections and an
% AND (exhaustive) stopping rule
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
v_b  = params(2); % v_b:  Rate in channel B ( multiplied by time scale )
p_ab = params(3); % p  :  Probability of sharing from a to b
p_ba = params(4); % p  :  Probability of sharing from b to a
c_a  = params(5); % A threshold
c_b  = params(6); % B threshold

%% Set up transition matrix

% P[1,1] = P(A Up & B Up)
p11 = v_a * v_b * (1-p_ab) * (1-p_ba);

% P[1,0] = P(A Up & B Same)
p10 = v_a * v_b * p_ab * (1-p_ba) +...
     v_a * (1-v_b) * (1-p_ab);
p10_b0 = v_a * (1-v_b) * p_ab; % added if b has no counts
p10_bc = v_a; % if b has reached criterion

% P[1,-1] = P(A Up & B Down)
p1_1 = v_a * (1-v_b) * p_ab;

% P[0,1] = P(A Same & B Up)
p01 = v_a * v_b * p_ba * (1-p_ab) +...
    (1-v_a) * v_b * (1-p_ba);
p01_a0 = (1-v_a) * v_b * p_ba; % added if a has no counts
p01_ac = v_b; % if a has reached criterion

% P[0,0] = P(A & B Same)
p00 = (1-v_a) * (1-v_b) +...
      v_a  * v_b * p_ab * p_ba;
p00_ac = 1-v_b; % if a has reached criterion 
p00_bc = 1-v_a; % if b has reached criterion

% P[-1,-1] = P(A Same & B Down)
p_1_1 = 0;

% P[-1,1] = P( A Down & B Up)
p_11 = (1-v_a) * v_b * p_ba;

% P[-1,0] = P(A Down & B Same)
p_10 = 0;

%%
T1_0  = diag([p00 * ones(c_b,1); p00_bc]) +...
         diag(p01 * ones(c_b,1), 1);              % Channel A at 0, B moves
T1_a0 = diag([p00 * ones(c_b,1);p00_bc]) +...
         diag((p01+p01_a0) * ones(c_b,1), 1);     % Channel A at intermediate counts, B moves
T1_ac = diag([p00_ac*ones(c_b,1);1]) +...
         diag(p01_ac*ones(c_b,1),1);              % Channel A has completed, B moves
    
T2_0  = diag([p1_1*ones(c_b-1,1);0], -1) + ...
        diag([p10+p10_b0;p10*ones(c_b-1,1);p10_bc]) + ...
        diag(p11 * ones(c_b,1), 1);
    
T3_0  = diag(p_11 * ones(c_b,1), 1);


T1 = mat2cell(repmat(T1_0, 1, c_a-1), c_b+1, repmat(c_b+1, 1, c_a-1));
T2 = mat2cell(repmat(T2_0, 1, c_a-1), c_b+1, repmat(c_b+1, 1, c_a-1));
% T3 = mat2cell(repmat(T3_0, 1, c_a-2), c_b+1, repmat(c_b+1, 1, c_a-2));
T3 = mat2cell(repmat(T3_0, 1, max(1,c_a-2)), c_b+1, repmat(c_b+1, 1, max(1,c_a-2)));

%% Build full transition matrix
T = blkdiag(T1_a0, T1{1:c_a-1}, T1_ac) + ...                        
    vertcat(blkdiag(horzcat(zeros(c_b+1,c_b+1), T2_0), T2{1:c_a-1}), ...
       zeros(c_b+1,(c_a+1)*(c_b+1))) + ...
    horzcat(blkdiag(vertcat(zeros(c_b+1,c_b+1), T3{1}), T3{1:c_a-2}, ...
       zeros(c_b+1,c_b+1)), zeros((c_b+1)*(c_a+1),c_b+1));

%% Internal transitions
k = size(T,1);             % Size of matrix
Q = T(1:(k-1), 1:(k-1));   % Transient states

%% Absorbing states
R = T(1:(k-1), k);         % Absorbing state vector

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