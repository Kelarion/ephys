function scatlab(varargin)
% Wrapper for scatter function which displays extra dimensions in the data
% cursor. 
% 
% scatlab(X,Y,Z,__)
%
% All extra arguments are passed directly into 'scatter' function.

fig = gcf;
% fig.DeleteFcn = 'doc datacursormode';
X = varargin{1};
Y = varargin{2};
lab = varargin{3};

scatter(X,Y,varargin{4:end})
dcm_obj = datacursormode(fig);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,lab})

%% datatip function
    function txt = myupdatefcn(~,event_obj,Z)
        % Customizes text of data tips
        pos = get(event_obj,'Position');
        
        ind=intersect(Find(X,pos(1),1e-10),Find(Y,pos(2),1e-10));
        if(length(ind)~=1)
            text='err';
        else
            text=num2str(Z(ind),4);
        end
        
        txt = {['X: ',num2str(pos(1))],...
            ['Y: ',num2str(pos(2))],...
            ['ID: ',text]};
    end

end

%% ----------------------------------------------------
function [ out ] = Find( vector, value ,precision)

if nargin < 3
    precision =   0.0001;
end

out=[];

for i=1:length(vector)
    
    if(abs(vector(i)-value)<precision)
        out=[out i];
    end

end

end