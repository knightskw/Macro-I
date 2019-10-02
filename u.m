% Shai Knight-Winnig 2016
% Helper file

function f = u(c)

global sigma CARA

if sigma == 1
    f = log(c);
elseif sigma == 'exp'
    f = -exp(-CARA*c)/CARA;
else
    f = c.^(1-sigma)/(1-sigma);
end
