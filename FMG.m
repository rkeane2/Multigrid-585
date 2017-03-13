function [u] = FMG(n,levels,probtype,mu1,mu2,w,relaxtype,resttype,iterV,outtype)
%note: only for probs 0,1,2,3,4
Fcycle(1).f = rhs(n,probtype,[]); 
Fcycle(1).n = n;
Fcycle(levels+1).f = [];
Fcycle(levels+1).n = [];

for i =1:levels
    Fcycle(i+1).n = n/(2^(i));
    Fcycle(i+1).f = interpop(Fcycle(i).n,Fcycle(i+1).n,resttype,probtype,Fcycle(i).f);
end

Fcycle(end).u = matrix(Fcycle(end).n,probtype)\Fcycle(end).f;

for i = levels:-1:1
    Fcycle(i).u = interpop(Fcycle(i).n,Fcycle(i+1).n,2,probtype,Fcycle(i+1).u);
    Fcycle(i).u = Vcycle(Fcycle(i).n,levels+1-i,Fcycle(i).u,probtype,mu1,mu2,w,relaxtype,resttype,iterV,0);
end

u = Fcycle(1).u;

if outtype ==1
if probtype ==1
    l_b = 1/2;   r_b = -1/2;
end
if probtype ==0
    l_b = 1/2;  r_b = 1;

end
if probtype ==2
    l_b = 1; r_b = 3;
end

if probtype == 0 || probtype ==1 || probtype ==2
u = [l_b;u;r_b];
end

if probtype ==3 
    u = reshape(u,(n-1),(n-1)); 
    u = [ones(1,n+1);ones(n-1,1),u,ones(n-1,1);ones(1,n+1)];
end
if probtype ==4
    u = v;
    u = reshape(u,(n-1),(n-1)); 
    u = [2*ones(1,n+1);2*ones(n-1,1),u,2*ones(n-1,1);2*ones(1,n+1)];
end
end


























