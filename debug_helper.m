
%just a code to help debugging
%can compare true with computed solution, 
%if using Vcycle with initial guess as true solution,
%it should stay within machine precision.
n = 256; 
probtype = 0;

A = matrix(n,probtype);
f = rhs(n,probtype);
% u_true = exact(n,probtype);
% u_true = u_true(2:end-1);
u_test = A\f;
 levels =3;
 w = 1.1; 
 relaxtype = 1;
 resttype = 1;
 u_test2 = Vcycle(n,levels,u_test,probtype,3,3,w,relaxtype,resttype,iters);
%norm(A*u_true-f)*sqrt(1/n)
%norm(A*u_test-f)*sqrt(1/n)
norm(u_test-u_test2(2:end-1))*sqrt(1/n)
norm(u_test-u_test2(2:end-1),inf)
%figure(1); clf;
plot([u_test,u_test2(2:end-1)])*sqrt(1/n);
%plot(u_true)

%% nonlinear debug

%testres = exact(64,4);
%testres = (testres(2:end-1,2:end-1));
%testres = zeros(63^2,1);
testres = testres(:);
testrhs = rhs(64,4,testres);
testjac = JM(testres);
testres = testres-testjac\testrhs;

%% nonlinear debug 2

u= u(2:end-1,2:end-1);
rhs(64,6,u(:))