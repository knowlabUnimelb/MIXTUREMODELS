function cdf = fac_conf(params, varargin)
% Poisson parallel interactive model with facilitatry connections and an
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
% shared, the other channel increments. If a channel generates an
% increment, you need to specify what happens to it (shared or not shared).
% Sharing can only occur when an increment is generated. 

% P[0,0] = P(A & B Same)
p00 = (1-v_a) * (1-v_d); % Neither channel increments
p00_ac = 1-v_d; % if a has reached criterion
p00_dc = 1-v_a; % if d has reached criterion

% P[1,0] = P(A Up & B Same)
p10 = v_a * (1-v_d) * (1-p_ad);  % A increments, D doesn't increment, A doesn't facilitate D
% p10_d0 = v_a * (1-v_d) * p_ad;  % added if d has no counts: A increments, D doesn't increment, A inhibits D (same as p1_1)
p10_dc = v_a; % if d has reached criterion

% P[0,1] = P(A Same & B Up)
p01 = (1-v_a) * v_d * (1-p_da);  % A doesn't increment, D increments, and D doesn't facilitate A
% p01_a0 = (1-v_a) * v_d * p_da; % added if a has no counts: A doesn't increment, D increments, D inhibits A
p01_ac = v_d; % if a has reached criterion 

% P[1,1] = P(A Up & B Up)
p11 = v_a * v_d * (1-p_ad) * (1-p_da) +... % A and D increment and nobody facilitates the other
      v_a * (1-v_d) * p_ad +...
      (1-v_a) * v_d * p_da;

% P[0,2] = P(A Same & B Up Twice)
p02 = 0;

% P[2,0] = P(A Up Twice & B Same)
p20 = 0;

% P[1,2] = P(A Up & B Up Twice)
p12 = v_a * v_d * p_ad * (1-p_da); % A increments, D increments, A facilitates D

% P[2,1] = P( A Up Twice & B Up)
p21 = v_a * v_d * p_da * (1-p_ad); % A increments, D increments, D facilitates A

% P[2,2] = P(A Up Twice & B Up Twice)
p22 = v_a * v_d * p_ad * p_da;


%%
T1_0  = diag([p00 * ones(c_d,1); p00_dc]) +...
         diag(p01 * ones(c_d,1), 1);              
T1_a0 = diag([p00 * ones(c_d,1); p00_dc]) +...
         diag(p01 * ones(c_d,1), 1);     
T1_ac = diag([p00_ac*ones(c_d,1);1]) +...
         diag(p01_ac*ones(c_d,1),1);         
        
T2_0  = diag([p10; p10*ones(c_d-1,1); p10_dc]) + ...
        diag([p11 * ones(c_d-1,1); p11+p12], 1) +...
        diag([p12*ones(c_d-1,1)], 2);
T2_ac = diag([p10; p10*ones(c_d-1,1); p10_dc]) +...
        diag([(p11+p21) * ones(c_d-1,1); p11+p12+p21+p22], 1) +...
        diag([(p12+p22)*ones(c_d-1,1)], 2);
    
    
T3_0  = diag([p21 * ones(c_d-1,1); p21+p22], 1) +...
        diag(p22 *ones(c_d-1,1), 2);

T1 = mat2cell(repmat(T1_0, 1, c_a-1), c_d+1, repmat(c_d+1, 1, c_a-1));
T2 = mat2cell(repmat(T2_0, 1, c_a-1), c_d+1, repmat(c_d+1, 1, c_a-1));
% T3 = mat2cell(repmat(T3_0, 1, c_a-2), c_b+1, repmat(c_b+1, 1, c_a-2));
T3 = mat2cell(repmat(T3_0, 1, max(1,c_a-2)), c_d+1, repmat(c_d+1, 1, max(1,c_a-2)));

%% Build full transition matrix
T = blkdiag(T1_a0, T1{1:c_a-1}, T1_ac) + ...                        
    vertcat(blkdiag(horzcat(zeros(c_d+1,c_d+1), T2_0), T2{1:c_a-2}, T2_ac), ...
       zeros(c_d+1,(c_a+1)*(c_d+1))) + ...
    vertcat(blkdiag(horzcat(zeros(c_a+1, 2*(c_d+1)), T3_0), T3{1:c_a-2}),...
       zeros(2*(c_a+1), (c_a+1)*(c_d+1)));
    
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