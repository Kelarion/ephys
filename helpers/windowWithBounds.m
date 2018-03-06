function inds = windowWithBounds(cntr,win,min,max)
% inds = windowWithBounds(cntr,win,min,max)
%
% makes a window, of a certain size, for indexing an array while respecting
% the limits of that array.
%
% cntr is the center index of the window
% win is of the form e.g. [-3 3]

istart = cntr + win(1);
ifin = cntr + win(2);
winSize = diff(win);

if istart >= min && ifin <= max
    inds = istart:ifin;
elseif istart <min
    inds = min:winSize;
else
    inds = (max-winSize):max;
end