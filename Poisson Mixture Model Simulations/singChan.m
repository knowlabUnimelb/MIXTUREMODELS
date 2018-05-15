function cdf = singChan(params, varargin)
% Poisson channel model
%
% The poisson distribution expresses the probability of a given number of 
%  events occurring in a fixed interval of time and/or space if these events 
%  occur with a known average rate and independently of the time since the 
%  last event
%
% This function implements a single channel of the parallel channel models 
% outlined in the Appendix of
% Johnson, Blaha, Houpt & Townsend (2010) JMP, 54, 53-72
%
% These single channels are useful for the computation of capacity

%% Optional arguments
optargs = {0:1000};
newVals = cellfun(@(x) ~isempty(x), varargin);% skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
[t] = optargs{:}; % Place optional args in memorable variable names

%% Parameters
v_a  = params(1);  %  v_a  :  accumulation rate 
c_a  = params(2);  %  c_a  : criterion for channel A

%% 
cdf = 1-poisscdf(c_a-1, v_a*t); % CDF
% S = poisscdf(c_a-1, v_a*t); % Survivor
