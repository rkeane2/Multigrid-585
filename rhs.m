function f = rhs(n,probtype)
%n = number of intervals
%probtype = model problem
%f = vector s.t. Au = f

if probtype ==1
    h = 1/n;
l_b = 1/2;
r_b = -1/2;
x = (1:(n-1))*h; x = x';
rhs_fun = @(x) 4*pi^2*sin(2*pi*x);
f = rhs_fun(x); 
f(1) = f(1)-l_b*(1/(h^2));
f(end) = f(end)-r_b*(1/(h^2));
end

if probtype ==0
    h = 1/n;
    l_b = 1/2;
r_b = 1;
x = (1:(n-1))*h; x = x';
rhs_fun = @(x) x.^2-(5/4)*x+3/8;
f = rhs_fun(x); 
f(1) = f(1)-l_b*(1/(h^2));
f(end) = f(end)-r_b*(1/(h^2));
end

if probtype ==2
    h = 1/n;
l_b = 1;
r_b = 3;
x = (1:(n-1))*h; x = x';
rhs_fun = @(x) -20+.5*120*pi*x.*cos(20*pi*x.^3)-.5*(60*pi*x.^2).^2.*sin(20*pi*x.^3);
f = rhs_fun(x); 
f(1) = f(1)-l_b*(1/(h^2));
f(end) = f(end)-r_b*(1/(h^2));
end

if probtype ==3
   h = 1/n;
   rhs_fun = @(x,y) x.^2+y.^2;
f = zeros((n-1)^2,1); 
for k=1:(n-1)^2
    count = 0;
    [i,j] = ind2sub([n-1,n-1],k);
    x=i*h ;
    y = j*h;
    f(k) = rhs_fun(x,y); 
    if i == 1 || i == n-1
        f(k) = f(k)-(1/h^2);
%         if count>0
%             f(k) = f(k) + (1/h^2);
%         end
%         count = count+1;
        
    end
    if j == 1 || j == n-1
        f(k) = f(k)-(1/h^2);
%         if count>0
%             f(k) = f(k) + (1/h^2);
%         end
    end
end
end




