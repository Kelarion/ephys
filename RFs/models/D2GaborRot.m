function F = D2GaborRot(x,xdata)
%% x = [Amp, x0, wx, y0, wy, wn, psi, fi]
% 
% Gaussian-modulated sine (Gabor), wn is in xdata units

xdatarot(:,:,1)= xdata(:,:,1)*cos(x(8)) - xdata(:,:,2)*sin(x(8));
xdatarot(:,:,2)= xdata(:,:,1)*sin(x(8)) + xdata(:,:,2)*cos(x(8));
x0rot = x(2)*cos(x(8)) - x(4)*sin(x(8));
y0rot = x(2)*sin(x(8)) + x(4)*cos(x(8));

f_renorm = x(6)*(2*pi);
s = sin(f_renorm*xdatarot(:,:,1) + x(7));

F = x(1)*s.*exp(   -(((xdatarot(:,:,1)-x0rot).^2)./(2*x(3)^2) + ((xdatarot(:,:,2)-y0rot).^2)./(2*x(5)^2) )    );