function F = D2_DoG_unequal_Rot(x,xdata)
%% x = [Amp1, x01, wx1, y01, wy1, Amp2, wx2, wy2, x02, y02, fi]
%
% Difference of two gaussians where means can be unequal

xdatarot(:,:,1)= xdata(:,:,1)*cos(x(11)) - xdata(:,:,2)*sin(x(11));
xdatarot(:,:,2)= xdata(:,:,1)*sin(x(11)) + xdata(:,:,2)*cos(x(11));

x01rot = x(2)*cos(x(11)) - x(4)*sin(x(11));
y01rot = x(2)*sin(x(11)) + x(4)*cos(x(11));
x02rot = x(9)*cos(x(11)) - x(10)*sin(x(11));
y02rot = x(9)*sin(x(11)) + x(10)*cos(x(11));

F = x(1)*exp(   -((xdatarot(:,:,1)-x01rot).^2/(2*x(3)^2) + (xdatarot(:,:,2)-y01rot).^2/(2*x(5)^2) )    ) ...
    - x(6)*exp(   -((xdatarot(:,:,1)-x02rot).^2/(2*x(7)^2) + (xdatarot(:,:,2)-y02rot).^2/(2*x(8)^2) )    );

