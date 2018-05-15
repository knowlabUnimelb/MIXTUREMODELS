function o = logit(x, type)

if nargin == 1
    type = 'standard';
end

switch type 
    case 'standard'
        o = log(x./(1 - x));
    case 'inverse'
        o = 1./(1 + exp(-x));
end