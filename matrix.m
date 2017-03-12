function A = matrix(n,probtype)
%n = number of intervals
%probtype = model problem
%A = matrix st Au = f
h = 1/n;
if probtype ==1 || probtype==0 ||probtype==2
main_diag = diag(2*ones(n-1,1));
super_diag = diag(-1*ones(n-2,1),1);
sub_diag = diag(-1*ones(n-2,1),-1);
A = main_diag+super_diag+sub_diag;
A = (-1/(h^2))*A; 
end

if probtype ==3
    A = zeros((n-1)^2,(n-1)^2);
    T = zeros(n-1,n-1); 
    I = eye(n-1); 
    T = T+-4*diag(ones(n-1,1))+diag(ones(n-2,1),1)+diag(ones(n-2,1),-1);
    for i = 1:n-1
A((n-1)*(i-1)+1:(n-1)*i,(n-1)*(i-1)+1:(n-1)*i) = T; 
    end
    for i=1:n-2
A((n-1)*(i-1)+1:(n-1)*i,(n-1)*(i)+1:(n-1)*(i+1)) = I; 
    end
for i=1:n-2
A((n-1)*(i)+1:(n-1)*(i+1),(n-1)*(i-1)+1:(n-1)*(i)) = I; 
end
A = A*(1/h^2);
end

