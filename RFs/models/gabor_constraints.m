function [c,ceq,gradc,gradceq] = gabor_constraints(x)
% [c,ceq,gradc,gradceq] = gabor_constraints(x)
%
% 

% x = [Amp, x0, wx, y0, wy, wn, psi, fi]'

c(1) = 1 - x(6)*x(3); % so, wn*wx >= 1, i.e. at least 1 cycle per 1std
ceq = [];

if nargout == 3
    gradc = [0;0;-x(6);0;0; -x(3);0;0]; % gradient of c(1)
    gradceq = [];
end


end