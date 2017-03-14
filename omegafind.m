function [opt_w] = omegafind(n,probtype,relaxtype)

% probtype = model problem to test
% relaxtype = 1 - full weighting, 0 - injection  
% n = 264; %number of grid points 

%make matrix M, look at eigenvalues. 

A = matrix(n,probtype); 
b = rhs(n,probtype); 
%v = linspace(1,3,n-1)'; 
if probtype == 0 || probtype ==1 || probtype ==2
v = zeros(n-1,1); 
end
if probtype ==3 || probtype ==4
    v = zeros((n-1)^2,1); 
end
if relaxtype ==0
    omega = .01:.01:1;
end
if relaxtype ==1
    omega = 1:.01:2;
end
 
%% method 1: look at eigenvalues of matrix 'M'. Not as good with multigrid. use method 2 for better multigrid
eigholder = zeros(length(A),length(omega));
for i = 1:length(omega)
    w = omega(i); 
[v,M] = WJac(A,b,[],w,0,relaxtype);
evals = eig(M); 
if relaxtype ==0 
evals = sort(abs(evals)); 
end
evals;
eigholder(:,i) = evals; 
end
%if relaxtype==0
eigholder = eigholder(1:round(end/2),:); 
%end
[~,w_index] = min(max(eigholder));
opt_w = omega(w_index); 

%note the optimal omega isn't consistent with numerical experiments

%% method 2: 
% testholder = zeros(length(omega),1); 
% for i = 1:length(omega)
%     w = omega(i);
%     test1 = WJac(A,b,v,w,5,relaxtype);
%     testholder(i) = norm(b-A*test1); 
% end
% 
% [~,w_index] = min(testholder); 
% opt_w = omega(w_index); 


