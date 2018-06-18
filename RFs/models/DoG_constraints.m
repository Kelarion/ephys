function [c,ceq,gradc,gradceq] = DoG_constraints(x)
% works for both DoG and DoG unequal, since A1 and A2 are found in the same
% place for those functions

% x = [A1,x01,wx1,y01,wy1,A2,wx2,wxy2,x02,y02,fi]'

c(1) = 0.0001 - x(1)*x(6); % amplitudes must have same sign (since A2 is subtracted)
% c(2) = x(5) - 5*x(3);
c(2) = x(8) - 5*x(7);
ceq = [];

if nargout == 3
    %     gradc = [-x(6) 0 0;
    %               0    0 0;
    %               0   -5 0;
    %               0    0 0;
    %               0    1 0;
    %              -x(1) 0 0;
    %               0    0 -5;
    %               0    0  1]; % gradient of c
    gradc = [-x(6) 0 ;
                  0     0;
                  0   0;
                  0    0;
                  0     0;
                 -x(1)  0;
                  0     -5;
                  0    1];
    gradceq = [];
end


end