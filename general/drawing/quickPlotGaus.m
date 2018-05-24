function quickPlotGaus(x,modfun)
% quickPlotGaus(x,modfun)
% 

faux = -130:2:130;
fauy = -33:2:33;
[X,Y] = meshgrid(fauy,faux);
xdata = Y;
xdata(:,:,2) = X;

F = modfun(x,xdata);
gcf;
imagesc(faux,-fauy,F')
