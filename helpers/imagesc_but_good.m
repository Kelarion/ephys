function imagesc_but_good(ax,x,y,mat,varargin)
% imagesc_but_good(ax,x,y,mat,varargin)
%
% Wrapper for imagesc that makes x and y actually specify the centers of
% each square, rather than only specifying the limits. But beware the
% xticks. Only use if the the axes dimensions are fixed.

imagesc(ax,mat,varargin{:})
set(ax,'XTick',1:size(mat,2));
set(ax,'YTick',1:size(mat,1));
ytix = yticks(ax); ytix = ytix(ytix == round(ytix));
xtix = xticks(ax); xtix = xtix(xtix == round(xtix));
set(ax,'XTickLabels',x(xtix))
set(ax,'YTickLabels',y(ytix))

end