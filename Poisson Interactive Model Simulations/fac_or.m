function cdf = fac_or(params, varargin)
% Poisson parallel interactive model with facilitatory connections and an
% OR (exhaustive) stopping rule
%
% The poisson distribution expresses the probability of a given number of 
%  events occurring in a fixed interval of time and/or space if these events 
%  occur with a known average rate and independently of the time since the 
%  last event
%
% This function implements the equation outlined in Box 2 of the Appendix
% of Johnson, Blaha, Houpt & Townsend (2010) JMP, 54, 53-72
%
% Note that: P(T1 < t OR T2 < t) = 1 - P(T1 > t AND T2 > t)

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
p = zeros(size(t)); % Preallocate

for x_a = 0:(c_a-1)
   pxa = poisspdf(x_a, t*v_a);
   for x_b = 0:(c_b-1)
      pxb = poisspdf(x_b, t*v_b);

      ba_facil = binocdf(min((c_a-x_a-1),x_b), x_b, p_ba);
      ab_facil = binocdf(min((c_b-x_b-1),x_a), x_a, p_ab);
      
      p = p + pxa .* pxb .* ab_facil .* ba_facil;
   end
end

cdf = 1 - p;