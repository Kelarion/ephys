function F = D2_DoG_Rot(x,xdata)
%% x = [Amp1, x0, wx1, y0, wy1, Amp2, wx2, wy2, fi]
% 
% Difference of two gaussians with equal means

xdatarot(:,:,1)= xdata(:,:,1)*cos(x(9)) - xdata(:,:,2)*sin(x(9));
xdatarot(:,:,2)= xdata(:,:,1)*sin(x(9)) + xdata(:,:,2)*cos(x(9));

x0rot = x(2)*cos(x(9)) - x(4)*sin(x(9));
y0rot = x(2)*sin(x(9)) + x(4)*cos(x(9));

F = x(1)*exp(   -((xdatarot(:,:,1)-x0rot).^2/(2*x(3)^2) + (xdatarot(:,:,2)-y0rot).^2/(2*x(5)^2) )    ) ...
    - x(6)*exp(   -((xdatarot(:,:,1)-x0rot).^2/(2*x(7)^2) + (xdatarot(:,:,2)-y0rot).^2/(2*x(8)^2) )    );