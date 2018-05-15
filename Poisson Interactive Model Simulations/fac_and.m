function cdf = fac_and(params, varargin)
% Poisson parallel interactive model with facilitatory connections and an
% AND (exhaustive) stopping rule
%
% The poisson distribution expresses the probability of a given number of 
%  events occurring in a fixed interval of time and/or space if these events 
%  occur with a known average rate and independently of the time since the 
%  last event
%
% This function implements the equation outlined in Box 1 of the Appendix
% of Johnson, Blaha, Houpt & Townsend (2010) JMP, 54, 53-72

%% Optional arguments
optargs = {0:1000};
newVals = cellfun(@(x) ~isempty(x), varargin);% skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
[t] = optargs{:}; % Place optional args in memorable variable names

%% Parameters
v_a  = params(1); %  v_a  :  accumulation rate in channel A
v_b  = params(2); %  v_b  :  accumulation rate  in channel B
p_ab = params(3); %  p    :  probability of interaction from channel A to B
p_ba = params(4); %  p    :  probability of interaction from channel B to A
c_a  = params(5); %  c_a  : criterion for channel A
c_b  = params(6); %  c_b  : criterion for channel B

%%
cdf = zeros(size(t)); % Preallocate cdf

%% Term 1
% Neither a nor b have enough to complete on their own but, rather, each
% requires some shared counts from the other channels
for x_a = 0:(c_a-1);
    for x_b = (c_a-x_a):(c_b-1) % This line eliminates the need to specify the min(x_b, c_a-x_a-1) in e.g., ba_facil
        
        % Within channel counts: total probability is the product of both channels 
        in_chan = poisspdf(x_a, t.*v_a) .*... % Probability of taking x_a steps in time t where the rate is t.*v_a
                  poisspdf(x_b, t.*v_b);      % Probability of taking x_b steps in time t where the rate is t.*v_b

        % Between channel interaction
        ba_facil = 1-binocdf(c_a-x_a-1, x_b, p_ba); % Here we only need to consider c_a-x_a-1 due to the construction of the x_b loop
        ab_facil = 1-binocdf(c_b-x_b-1, x_a, p_ab);

        cdf = cdf + in_chan .* ba_facil .* ab_facil;
    end
end

%% Term 2
% B completes on its own but A needs help
x_b_complete = 1-poisscdf(c_b-1,t.*v_b); % Channel B completes through the poisson channel along
for x_a = 0:(c_a-1) % Channel A requires shared counts from B to finish
    % Within channel counts
    in_chan = poisspdf(x_a,t.*v_a);

    % Between channel interaction
    ba_facil = 1-binocdf(c_a-x_a-1, c_b, p_ba);

    cdf = cdf + x_b_complete .* in_chan .* ba_facil;
end

%% Term 3
% A completes on its own but B needs help
x_a_complete = 1-poisscdf(c_a-1,t.*v_a);
for x_b = 0:(c_b-1)
     % Within channel counts
     in_chan = poisspdf(x_b, t.*v_b);

     % Between channel interaction
     ab_facil = 1-binocdf(c_b-x_b-1, c_a, p_ab);

     cdf = cdf + x_a_complete .* in_chan .* ab_facil;
end

% Both channels complete on their own
cdf = cdf + x_b_complete .* x_a_complete;
% Same as: cdf = cdf + (1-poisscdf(c_b-1,t.*v_b)) .* (1-poisscdf(c_a-1,t.*v_a) );