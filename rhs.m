function f = rhs(n,probtype,v)
%n = number of intervals
%probtype = model problem
%f = vector s.t. Au = f
%v = current iter for probtype = 6
%if probtype =6 it returns the current residual = f-Av

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
   rhs_fun = @(x,y) (x-.5).^2+y.^2;
f = zeros((n-1)^2,1); 
for k=1:(n-1)^2
    count = 0;
    [j,i] = ind2sub([n-1,n-1],k);
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





if probtype ==4 %experiment with your favorite function :) 
   h = 1/n;
   rhs_fun = @(x,y) 120*pi*x.*(sin(pi*(x+1).^3))-60*pi*y.*cos(pi^3.*(4*y-2).^2);
%rhs_fun = @(x,y) 0;
   f = zeros((n-1)^2,1); 
for k=1:(n-1)^2
    [j,i] = ind2sub([n-1,n-1],k);
    x=i*h ;
    y = j*h;
    f(k) = rhs_fun(x,y); 
    if i == 1 %top
        f(k) = f(k)-(2/h^2);       
    end
    if j == 1 %left
        f(k) = f(k)-(2/h^2);
    end
    if i == n-1 %bottom
        f(k) = f(k)-(2/h^2);       
    end
    if j == n-1 %right
        f(k) = f(k)-(2/h^2);
    end
end
end

if probtype ==5
   h = 1/n;
   rhs_fun = @(x,y) (x-.5).^2+1000*sin(48*pi*y.^2);
f = zeros((n-1)^2,1); 
for k=1:(n-1)^2
    count = 0;
    [j,i] = ind2sub([n-1,n-1],k);
    x=i*h ;
    y = j*h;
    f(k) = rhs_fun(x,y); 
    if i == 1 %top
        f(k) = f(k)-(1/h^2);       
    end
    if j == 1 %left
        f(k) = f(k)-(1/h^2);
    end
    if i == n-1 %bottom
        f(k) = f(k)-(1/h^2);       
    end
    if j == n-1 %right
        f(k) = f(k)-(1/h^2);
    end
end
end

if probtype ==6
    gam = 10; %can change
    h = 1/n;
    x = (1:n-1)*h;
    [X,Y] = meshgrid(x,x); 
    rhs_fun = @(x,y) 2*((x-x.^2)+(y-y.^2))+gam*(x-x.^2).*(y-y.^2).*exp((x-x.^2).*(y-y.^2));
    
    f = rhs_fun(X,Y);
    f = f(:); 
    v = reshape(v,n-1,n-1);
    v = [zeros(1,n+1);zeros(n-1,1),v,zeros(n-1,1);zeros(1,n+1)];
    A = zeros(n+1,n+1);
    for i=2:n
        for j=2:n
            A(i,j) = (1/(h^2))*(4*v(i,j)-v(i-1,j)-v(i+1,j)-v(i,j-1)-v(i,j+1))+...
                gam*v(i,j)*exp(v(i,j));
        end
    end
    A = A(2:end-1,2:end-1); A = A(:);
    f = f-A; %f = f-A;
end




