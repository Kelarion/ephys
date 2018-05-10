function t = nohomopois(intens,T)
% t = nohomopois(intens,T)
%
% example of generating a nonhomogeneous poisson process 
% on [0,T] with intensity function intens
%
% e.g. t = nonhomopp(@(x)(50*sin(x)),10*pi);
% 
% Made by some guy

if isvector(T)
    lambda = interp1(T,intens(T),1:0.001:max(T));
    intens = @(t,lambda)(lambda(findNearestPoint(t,1:0.001:T)));
    T = max(T);
else
    lambda = 1;
    intens = @(t,lambda)(lambda*intens(t));
end

x = 0:.01:T;
m = intens(x,lambda);
m2 = max(m); % generate homogeneouos poisson process
u = rand(1,ceil(1.5*T*m2));
t = cumsum(-(1/m2)*log(u)); %points of homogeneous pp
t = t(t<T); n=length(t); % select those points less than T
m = intens(t,lambda); % evaluates intensity function
t = t(rand(1,n)<m/m2); % filter out some points