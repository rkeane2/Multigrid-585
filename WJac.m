function [v,M] = WJac(A,f,v,w,iter,relaxtype)
%input:
%v = initial guess
if relaxtype ==0
D = diag(diag(A));
R = A-D; 
c = w*(D\f);
M = -w*(D\R)+(1-w)*eye(length(A));

for i=1:iter
    v = M*v+c;
end
end

if relaxtype ==1
    D = diag(diag(A));
L = tril(A)-D;
U = triu(A)-D; 
M = -(D+w*L)\(w*U+(w-1)*D);
    c = (D+w*L)\(w*f);

for i=1:iter
    v = M*v+c;
end
end
