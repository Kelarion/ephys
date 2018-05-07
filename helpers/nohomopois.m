function t = nohomopois(intens,T)
% t = nohomopois(intens,T)
%
% example of generating a nonhomogeneous poisson process 
% on [0,T] with intensity function intens
%
% e.g. t = nonhomopp(@(x)(50*sin(x)),10*pi);
% 
% Made by some guy

x = 0:.0001:T;
m = intens(x);
m2 = max(m); % generate homogeneouos poisson process
u = rand(1,ceil(1.5*T*m2));
t = cumsum(-(1/m2)*log(u)); %points of homogeneous pp
t = t(t<T); n=length(t); % select those points less than T
m = intens(t); % evaluates intensity function
t = t(rand(1,n)<m/m2); % filter out some pointst