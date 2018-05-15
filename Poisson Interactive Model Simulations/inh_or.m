function cdf = inh_or(params, varargin)
% Poisson parallel interactive model with inhibitory connections and an
% AND (exhaustive) stopping rule
%
% The poisson distribution expresses the probability of a given number of 
%  events occurring in a fixed interval of time and/or space if these events 
%  occur with a known average rate and independently of the time since the 
%  last event
%
% This function implements the transition matrix method outlined in Table A.2 
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
% P( A & B Same)
p00 = v_a * v_b * p_ab * p_ba +...
      (1-v_a) * (1-v_b);

% P( A Up & B Same)
p10 = v_a * v_b * p_ab * (1-p_ba) +...
      v_a * (1-v_b) * (1-p_ab);
% p10_b0 = v_a * (1-v_b) * p_ab; % added if b has no counts
  
% P( A Same & B Up)
p01 = v_a * v_b * p_ba * (1-p_ab) +...
      (1-v_a) * v_b * (1-p_ba);
% p01_a0 = (1-v_a) * v_b * p_ba; % added if a has no counts

% P( A Up & B Down)
p1_1 = v_a * (1-v_b) * p_ab;

% P( A Down & B Up)
p_11 = (1-v_a) * v_b * p_ba;

% P( A Up & B Up)
p11 = v_a * v_b * (1-p_ab) * (1-p_ba);

%%
T1_a0 = diag([p00 * ones(c_b,1); 1]) +...
        diag((p01+p_11) * ones(c_b,1), 1);
T1_0  = diag([p00 * ones(c_b,1); 1]) +...
         diag(p01 * ones(c_b,1), 1);
T1_ac = eye(c_b+1);

T2_0 = diag([p1_1*ones(c_b-1,1);0], -1) + ...
        diag([p10+p1_1; p10*ones(c_b-1,1); 0]) + ...
        diag(p11 * ones(c_b,1), 1);

% T2_0  = diag([p1_1*ones(c_b-1,1);0], -1) + ...
%         diag([p10+p10_b0;p10*ones(c_b-1,1);0]) + ...
%         diag(p11 * ones(c_b,1), 1);
    
T3_0  = diag(p_11 * ones(c_b,1), 1);


T1 = mat2cell(repmat(T1_0, 1, c_a-1), c_b+1, repmat(c_b+1, 1, c_a-1));
T2 = mat2cell(repmat(T2_0, 1, c_a-1), c_b+1, repmat(c_b+1, 1, c_a-1));
T3 = mat2cell(repmat(T3_0, 1, c_a-2), c_b+1, repmat(c_b+1, 1, c_a-2));

%% Build full transition matrix
T = blkdiag(T1_a0, T1{1:c_a-1}, T1_ac) + ...
    vertcat(blkdiag(horzcat(zeros(c_b+1,c_b+1), T2_0), T2{1:c_a-1}), ...
       zeros(c_b+1,(c_a+1)*(c_b+1))) + ...
    horzcat(blkdiag(vertcat(zeros(c_b+1,c_b+1), T3{1}), T3{1:c_a-2}, ...
       zeros(c_b+1,c_b+1)), zeros((c_b+1)*(c_a+1),c_b+1));

%% Internal transitions
% Absorbing states
k = size(T,1);                      % Size of matrix
indR = [(c_b+1)*(1:c_a) (k-c_b):k]; % Index of absorbing states
indQ = setdiff(1:k, indR);          % Index of transient states
Q= T(indQ, indQ);                   % Transient states
R= T(indQ, indR);                   % Absorbing states

%%
k = size(Q,1);
cdf = zeros(size(t));
I = eye(k);
for i = 1:size(t,2)
    n = t(i);
    C_N = sum(inv(I - Q) * (I - Q ^(n+1)) * R, 2);
    cdf(i) = C_N(1);
end