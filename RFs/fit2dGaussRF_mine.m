

function [x, fitResp] = fit2dGaussRF_mine(xcoords, ycoords, resp, modfun,params)
% resp = your measure of the responsiveness of the neuron at each point in
% space
% x = [amplitude, xCenter, xStdev, yCenter, yStdev, rotation, baseline]
%
% specify what model you want, and if you do then also supply params, which
% is [x0;lb;ub]

    use2dgWithBsl = 0;
    opts = optimset('Display','off');
    
    if nargin<4, modfun = @D2GaussFunctionRot; end
    supplied = ~(nargin<5);
    
    [x, y] = meshgrid(xcoords, ycoords);
    xdata = zeros(size(x,1), size(x,2), 2);
    xdata(:,:,1) = x;
    xdata(:,:,2) = y;
    
    % x0 is the initial guess
    if supplied
        x0 = params(1,:);
    else
        [maxY, maxX] = find(resp==max(resp(:)),1);
        x0 = [1,maxX,mean(diff(xcoords))*5,maxY,mean(diff(ycoords))*5,0];
    end
    
    if use2dgWithBsl
        
        x0(7) = 0;

        lb = [0,min(xcoords),0,min(ycoords),0,-pi/4, 0];
        ub = [realmax('double'),max(xcoords),(max(xcoords))^2,max(ycoords),(max(ycoords))^2,pi/4, max(resp(:))];
        [x,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunctionRotWithBsl,x0,xdata,resp,lb,ub,opts);

        fitResp = D2GaussFunctionRotWithBsl(x, xdata);
        
    else
        
        if supplied
            lb = params(2,:);
            ub = params(3,:);
        else
            lb = [0,min(xcoords),0,min(ycoords),0,-pi/4];
            ub = [realmax('double'),max(xcoords),(max(xcoords))^2,max(ycoords),(max(ycoords))^2,pi/4];
        end
        [x,resnorm,residual,exitflag] = lsqcurvefit(modfun,x0,xdata,resp,lb,ub,opts);

        fitResp = modfun(x, xdata);
        
    end
    
%     if makeplots
%         figure;
%         subplot(3,1,1);
%         imagesc(xcoords, ycoords, resp);
%         colormap hot;
%         
%         subplot(3,1,2);
%         imagesc(xcoords, ycoords, fitResp);
%         
%         subplot(3,1,3);
%         contour(fitResp,1,'Color', 'k', 'LineWidth', 2.0);
%         set(gca, 'YDir', 'reverse');
%     end
end