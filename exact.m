function [u,x] = exact(n,probtype)
%n = number of intervals
%probtype = model problem
%u = exact solution

if probtype ==1
    h = 1/n; 
x = (1:(n-1))*h; x = x';x = [0;x;1];
exact_fun = @(x) -sin(2*pi*x)-x+1/2;
u = exact_fun(x); 
end
if probtype ==0
    h = 1/n; 
    x = (1:(n-1))*h; x = x';x = [0;x;1];
exact_fun = @(x) (1/12)*x.^4-(5/24)*x.^3+(3/16)*x.^2+(21/48)*x+1/2;
u = exact_fun(x); 
end
if probtype ==2
    h = 1/n; 
    x = (1:(n-1))*h; x = x';x = [0;x;1];
exact_fun = @(x) 1+12*x-10*x.^2+.5*sin(20*pi*x.^3);
u = exact_fun(x); 
end

if probtype ==3
    tar_dir = 'C:\Users\ronan\Documents\matlab code\mathematica code\UW etc\WI 2017\585 hw4 hw5 hw6\chebfun-master\chebfun-master';
    parent_dir = 'C:\Users\ronan\Documents\matlab code\mathematica code\UW etc\WI 2017\585 hw4 hw5 hw6\585 final';
    h = 1/n;
    x = (0:n)*h; 
    [X,Y] = meshgrid(x,x); 
    x = [X,Y];
    cd(tar_dir); 
    rhs_fun = chebfun2(@(x,y)x.^2+y.^2,[0 1 0 1]);
    A = chebop2(@(u) diff(u,2,1)+diff(u,2,2),[0 1 0 1]);
    A.lbc = 1; A.rbc = 1; A.ubc = 1; A.dbc = 1;
    u_cheb = A\rhs_fun;
    u = u_cheb(X,Y); 
    cd(parent_dir);
end

if probtype ==4
    tar_dir = 'C:\Users\ronan\Documents\matlab code\mathematica code\UW etc\WI 2017\585 hw4 hw5 hw6\chebfun-master\chebfun-master';
    parent_dir = 'C:\Users\ronan\Documents\matlab code\mathematica code\UW etc\WI 2017\585 hw4 hw5 hw6\585 final';
    h = 1/n;
    x = (0:n)*h; 
    [X,Y] = meshgrid(x,x); 
    x = [X,Y];
    cd(tar_dir); 
    rhs_fun = chebfun2(@(x,y)120*pi*x.*(sin(pi*(x+1).^3))-60*pi*y.*cos(pi^3.*(4*y-2).^2),[0 1 0 1]);
    A = chebop2(@(u) diff(u,2,1)+diff(u,2,2),[0 1 0 1]);
    A.lbc = 2; A.rbc = 2; A.ubc = 2; A.dbc = 2;
    u_cheb = A\rhs_fun;
    u = u_cheb(X,Y); 
    cd(parent_dir);
end
    

if probtype ==5
    h = 1/n;
    x = (0:n)*h; 
    [X,Y] = meshgrid(x,x);
    exact_fun = @(x,y) (x-x.^2).*(y-y.^2);
    x = [X,Y];
    u = exact_fun(X,Y);
end
    
    
