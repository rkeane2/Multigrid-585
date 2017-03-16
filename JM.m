function Jmat = JM(v)
%jacobian for nonlinear problem
gam = 10; %can change but remember to change in rhs as well
m = length(v); %m = (n-1)^2
Jmat = zeros(m,m); 
n = sqrt(m); %n = n-1
h = 1/(n+1); 

B = -eye(n,n)*(1/(h^2));
superdiag = diag(-ones(m-1,1),1)*(1/(h^2));
subdiag = diag(-ones(m-1,1),-1)*(1/(h^2));
diag_f =@(x) 4/(h^2)+gam.*x.*exp(x);
maindiag = diag(diag_f(v));
Jmat = Jmat+maindiag+subdiag+superdiag;
for i=1:n-1
Jmat((n)*(i-1)+1:(n)*i,(n)*(i)+1:(n)*(i+1)) = B; 
    end
for i=1:n-1
Jmat((n)*(i)+1:(n)*(i+1),(n)*(i-1)+1:(n)*(i)) = B; 
end






























