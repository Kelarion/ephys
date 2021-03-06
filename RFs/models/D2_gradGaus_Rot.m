function F = D2_gradGaus_Rot(x,xdata)
%% x = [Amp, x0, wx, y0, wy, fi]
%
% Vector gradient of a gaussian

xdatarot(:,:,1)= xdata(:,:,1)*cos(x(6)) - xdata(:,:,2)*sin(x(6));
xdatarot(:,:,2)= xdata(:,:,1)*sin(x(6)) + xdata(:,:,2)*cos(x(6));

x0rot = x(2)*cos(x(6)) - x(4)*sin(x(6));
y0rot = x(2)*sin(x(6)) + x(4)*cos(x(6));

a = -( (xdatarot(:,:,1)-x0rot)*x(5)^2 - (xdatarot(:,:,2)-y0rot)*x(3)^2)/((x(3)^2)*x(5)^2);

F = a.*(x(1)*exp(   -((xdatarot(:,:,1)-x0rot).^2/(2*x(3)^2) + (xdatarot(:,:,2)-y0rot).^2/(2*x(5)^2) )    ) );
