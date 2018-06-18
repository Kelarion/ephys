function F = D2_LoG_Rot(x,xdata)
%% x = [Amp, x0, wx, y0, wy, fi]
%
% Laplacian of a gaussian

xdatarot(:,:,1)= xdata(:,:,1)*cos(x(6)) - xdata(:,:,2)*sin(x(6));
xdatarot(:,:,2)= xdata(:,:,1)*sin(x(6)) + xdata(:,:,2)*cos(x(6));

x0rot = x(2)*cos(x(6)) - x(4)*sin(x(6));
y0rot = x(2)*sin(x(6)) + x(4)*cos(x(6));

% this is super gross but it is what it is
ax = ( x0rot*(2*xdatarot(:,:,1)-x0rot) - (x(3)^2 - xdatarot(:,:,1).^2))/(x(3)^4);
ay = (y0rot*(2*xdatarot(:,:,2)-y0rot) - (x(5)^2 - xdatarot(:,:,2).^2))/ (x(5)^4);

F = (ax+ay).*(x(1)*exp(   -((xdatarot(:,:,1)-x0rot).^2/(2*x(3)^2) + (xdatarot(:,:,2)-y0rot).^2/(2*x(5)^2) )    ) );

