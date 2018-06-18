function bounds = setBounds(mdl,rf,xPos,yPos,init_std)


[~, maxX] = find(rf==max(rf(:)),1);
[~, minX] = find(rf==min(rf(:)),1);
foo = [maxX minX]; % choose hemifield based on strongest peak
bar = [abs(max(rf(:))) abs(min(rf(:)))];
maxX = foo(bar == max(bar));

mid = mean(xPos);
hmfld = sign(xPos(maxX) - mid) < 0; % true if left hemifield

if hmfld % restrict search to hemifield that maximum is in
    rf = rf(:,xPos<mid);
    xPos = xPos(xPos<mid);
    [maxY, maxX] = find(rf==max(rf(:)),1);
    [minY, minX] = find(rf==min(rf(:)),1);
    x_min = min(xPos); % bounds
    x_max = mid;
else
    rf = rf(:,xPos>mid);
    xPos = xPos(xPos>mid);
    [maxY, maxX] = find(rf==max(rf(:)),1);
    [minY, minX] = find(rf==min(rf(:)),1);
    x_min = mid; % bounds
    x_max = max(xPos);
end
if mdl.fitIpsi 
    x_min = x_min + (~hmfld)*20; % shift minimum right if in right hemifield
    x_max = x_max - (hmfld)*20; % shift maximum left if in left hemifield
end

maxZ = (max(rf(:))-mean(rf(:)))./std(rf(:));
minZ = (min(rf(:))-mean(rf(:)))./std(rf(:));
if abs(minZ)>maxZ
    % if peak is negative, guess is at minimum value
    mu_x1 = xPos(minX);
    mu_y1 = yPos(minY);
else % otherwise the maximum
    mu_x1 = xPos(maxX);
    mu_y1 = yPos(maxY);
end

y_min = min(yPos);
y_max = max(yPos);
sx_max = range(xPos)/2;
sx_min = abs(min(diff(xPos)));
sy_max = range(yPos)/2;
sy_min = abs(min(diff(yPos)));

w_max = 0.3/abs(min(diff(xPos))); % avoid ailiasing

maxamp = 2*max(rf(:));
minamp = 2*min(rf(:));

std_guess = abs(mean(diff(xPos)))*init_std;

switch lower(mdl.String)
    case 'gaussian'
        lb = [minamp, x_min, sx_min, y_min, sy_min, -pi/4];
        ub = [maxamp,x_max,sx_max, y_max, sy_max,pi/4];
        
        % x = [A1,x0,wx,y0,wy,fi]'
        % 0.2 < wx1/wxy1 < 5 (can't be too oblong)
        A_ineq = full(sparse([3 3,5 5],[3 5,3 5],[1 -5,-5 1],6,6));
        b_ineq = zeros(6,1);
        nlin = [];
        
        x0 = [1,mu_x1,std_guess,mu_y1,std_guess,0];
    case 'dog'
        lb = [minamp, x_min, sx_min, y_min, sy_min, -maxamp, sx_min, sy_min, -pi];
        ub = [maxamp, x_max, sx_max, y_max, sy_max, ...
            -minamp, 0.8*sx_max, 0.8*sy_max, pi];
        
        % x = [A1,x01,wx1,y01,wy1,A2,wx2,wxy2,fi]'
        % constraints are: wx1 >= wx2, wx1 >= wy1,  
        % wy1 >= wy2, wx2 >= wy2 and A2>= A1
        %         A_ineq = full(sparse([3 3,5 5,7 7,8 8, 1 1],[3 5,3 5,3 7,5 8, 1 6],[-5 1,5 -1,-1 1,-1 1,1 -1],9,9));
        rws = [3 3, 5 5, 7 7, 8 8, 1 1];
        cls = [3 7, 3 5, 7 8, 5 8, 1 6];
        vls = [-1 1,-1 1,-1 1,-1 1,1 -1];
        A_ineq = full(sparse(rws,cls,vls,9,9));
        b_ineq = zeros(9,1); % A*x <= b are the constraints
        % additionally, the amplitudes should be opposite sign
        nlin = @DoG_constraints;
        
        x0 = [1, mu_x1, std_guess, mu_y1, std_guess, ...
            1.4, 0.6*std_guess, 0.6*std_guess, 0];
    case 'doug'
        lb = [minamp, x_min, sx_min, y_min, sy_min, ...
            -maxamp, sx_min,sy_min,x_min,y_min -pi/4];
        ub = [maxamp, x_max, sx_max, y_max, sy_max, ...
            -minamp, sx_max, sy_max, x_max, y_max, pi/4];
        
        mu_x1 = mean(xPos([maxX minX])); mu_x2 = mu_x1;
        mu_y1 = mean(yPos([maxY minY])); mu_y2 = mu_y1;
        
        % x = [A1,x01,wx1,y01,wy1,A2,wx2,wxy2,x02,y02,fi]'
        % we want wx1 > wx2 and wy1 > wy2, A2 >= A1, and
        % 0.2 < wx1/wxy1 < 5 (which constrains wx2/wy2 also)
        rws = [2 2 2,   3 3, 4 4 4,  5 5, 7 7, 8 8, 9 9 9, 10 10 10, 1 1];
        cls = [2 3 9,   3 7, 4 5 10, 3 5, 7 8, 5 8, 2 3 9, 4 5 10, 1 6];
        vls = [-1 -1 1,-1 1,-1 -1 1,-1 1,-1 1,-1 1,1 -1 -1,1 -1 -1,1 -1];
        A_ineq = full(sparse(rws,cls,vls,11,11));
        b_ineq = zeros(11,1); % A*x <= b are the constraints
        % additionally, the amplitudes should be opposite
        nlin = @DoG_constraints;
        
        x0 = [1, mu_x1, std_guess, mu_y1, std_guess, ...
            1, 0.6*std_guess, 0.6*std_guess, mu_x2, mu_y2, 0];
    case 'gabor'
        mu_x1 = xPos(maxX);
        mu_y1 = yPos(maxY);
        fi_guess = -atan((yPos(maxY) - yPos(minY))/(xPos(maxX) - xPos(minX)));
        
        lb = [minamp, x_min, 0, y_min, 0, 0, -pi/2, -pi];
        ub = [maxamp,x_max,sx_max, y_max, sy_max, ...
            w_max, pi/2, pi];
        
        % x = [A1,x0,wx,y0,wy,wn,psi,fi]'
        % 1/3 < wx1/wxy1 < 3 (can't be too oblong)
        A_ineq = full(sparse([3 3,5 5],[3 5,3 5],[1 -3,-3 1],8,8));
        b_ineq = zeros(8,1);
        % and sinusoid needs at least 1 cycle per std of mask
        nlin = @gabor_constraints;
        
        x0 = [1, mu_x1, 1.3*std_guess, mu_y1, std_guess, ...
            0.1/abs(min(diff(xPos))), 0, fi_guess];
    case 'flat'
        lb = [];
        ub = [];
        x0 = [];
        A_ineq = [];
        b_ineq = [];
        nlin = [];
    case 'rf'
        lb = [];
        ub = [];
        x0 = [];
        A_ineq = [];
        b_ineq = [];
        nlin = [];
end


bounds.x0 = x0;
bounds.lb = lb;
bounds.ub = ub;
bounds.A_ineq = A_ineq;
bounds.b_ineq = b_ineq;
bounds.nlin = nlin;

end

