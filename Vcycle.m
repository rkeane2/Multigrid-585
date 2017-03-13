function [u] = Vcycle(n,levels,v,probtype,mu1,mu2,w,relaxtype,resttype,iters,outtype)

% n = number of intervals on finest grid

%levels = number of levels to use on V cycle. 
%each level has n/2 intervals where n is the number of 
%intervals for the next finer level

%v is the initial guess

%probtype is the type of problem to solve
%probtype = 0: u'' = x^2-(5/4)x+3/8; u(0) = 1/2; u(1) = 1;
%probtype = 1: u'' = 4*pi^2*sin(2*pi*x); u(0) = 1/2; u(1) = 1/2;

%mu1  %number of times to relax on way down
%mu2  %number of times to relax on way up

%w value of omega for Jacobi/G-S. w =1 corresponds to regular Jacobi/G-S

%relax type is method to use for relaxation
%relaxtype = 0: weighted Jacobi
%relaxtype = 1: weighted G-S (SOR)

%resttype = type of restriction operation:
%resttype = 0: injection 
%resttype = 1: full weighting

%outtype = 0; output for FMG
%outtype = 1; output for plotting


if isempty(v)
    v = zeros(n-1,1); %default initial guess to zeros
end
if isempty(mu1)
    mu1 = 3; %default 3 iterations of relaxation
end
if isempty(mu2)
    mu2 = 3; %default 3 iterations of relaxation
end
if isempty(w)
    w = 1; %default w = 1
end
if isempty(relaxtype)
    relaxtype = 1; %default to G-S
end
if isempty(resttype)
    resttype = 1; %default to full weighting
end
if isempty(iters)
    iters = 1; %default to 1 iteration
end

Vgrid(levels+1).n = [];
Vgrid(levels+1).A = [];
Vgrid(levels+1).f = [];
Vgrid(1).f = rhs(n,probtype,v); 
for i =1:levels+1
Vgrid(i).n = n/(2^(i-1));
if probtype <6
Vgrid(i).A = matrix(Vgrid(i).n,probtype);
continue
end
if probtype ==6
    if i==1
        Vgrid(i).v = v;
    end
    Vgrid(i).A = JM(Vgrid(i).v);
    Vgrid(i+1).v = interpop(Vgrid(i).n,(Vgrid(i).n)/2,resttype,probtype,Vgrid(i).v);
end
end
if probtype ==1
    l_b = 1/2;
    r_b = -1/2;
%     h=1/n;
%     x = (1:(n-1))*h; x = x';x = [0;x;1];
end
if probtype ==0
    l_b = 1/2;
    r_b = 1;
%     h = 1/n;
%     x = (1:(n-1))*h; x = x';x = [0;x;1];
end
if probtype ==2
    l_b = 1;
    r_b = 3;
%     h = 1/n;
%     x = (1:(n-1))*h; x = x';x = [0;x;1];
end

% if probtype ==3
% %     h = 1/n;
% %     x = (0:n)*h; 
% %     [X,Y] = meshgrid(x,x); 
% %     x = [X,Y]; %grid for 2-d concatenated
% end

for k=1:iters
    
Vgrid(1).u = WJac(Vgrid(1).A,Vgrid(1).f,v,w,mu1,relaxtype); %relax initial guess mu1 times
Vgrid(2).f = Vgrid(1).f-Vgrid(1).A*Vgrid(1).u; %compute residual

%h = 1/n; norm(Vgrid(2).f)*sqrt(h^2) %print out residuals at each grid

for i = 1:levels
%     if probtype<3
% Vgrid(i+1).f = interpop(Vgrid(i).n,Vgrid(i+1).n,resttype,probtype,Vgrid(i+1).f)*Vgrid(i+1).f;  %compute coarse residual
%     end
%     if probtype>2
        Vgrid(i+1).f = interpop(Vgrid(i).n,Vgrid(i+1).n,resttype,probtype,Vgrid(i+1).f); %coarse residual
%     end
if i==levels
    Vgrid(i+1).u=Vgrid(i+1).A\Vgrid(i+1).f;
    break
end

    if probtype >2
    init_guess = zeros((Vgrid(i+1).n-1)^2,1);
else
    init_guess = zeros(Vgrid(i+1).n-1,1); 
end

Vgrid(i+1).u = WJac(Vgrid(i+1).A,Vgrid(i+1).f,init_guess,w,mu1,relaxtype); %compute error = "sol"
Vgrid(i+2).f = Vgrid(i+1).f-Vgrid(i+1).A*Vgrid(i+1).u; %compute new residual

%norm(Vgrid(i+2).f)*sqrt(h^2) %print out residuals at each grid

end

for i=levels:-1:1
%     if probtype<3
%     Vgrid(i).u = Vgrid(i).u+interpop(Vgrid(i).n,Vgrid(i+1).n,2,probtype,Vgrid(i+1).u)*Vgrid(i+1).u;
%     end
%     if probtype>2
    Vgrid(i).u = Vgrid(i).u+interpop(Vgrid(i).n,Vgrid(i+1).n,2,probtype,Vgrid(i+1).u);
%     end
    Vgrid(i).u = WJac(Vgrid(i).A,Vgrid(i).f,Vgrid(i).u,w,mu2,relaxtype);
end

v = Vgrid(1).u;
if outtype ==0
    u = v;
end
end

if outtype ==1
if probtype == 0 || probtype ==1 || probtype ==2
u = [l_b;v;r_b];
end

if probtype ==3
    u = v;
    u = reshape(u,(n-1),(n-1)); 
    u = [ones(1,n+1);ones(n-1,1),u,ones(n-1,1);ones(1,n+1)];
end
if probtype ==4
    u = v;
    u = reshape(u,(n-1),(n-1)); 
    u = [2*ones(1,n+1);2*ones(n-1,1),u,2*ones(n-1,1);2*ones(1,n+1)];
end
if probtype ==5
    u = v; h = 1/n;
    u = reshape(u,(n-1),(n-1)); 
    u = [1.*ones(1,n+1);1*ones(n-1,1),u,1*ones(n-1,1);1.*ones(1,n+1)];
end
end









