
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