%plot the error from SOR with different omega on 
probtype = 2;
relaxtype =1;
n = 264; 
h = 1/n;
iters = 20; 
A = matrix(n,probtype); 
f = rhs(n,probtype); 
u_true = exact(n,probtype); u_true = u_true(2:end-1); 
v = zeros(n-1,1); 
ploterror = zeros(iters+1,4); 
ploterror(1,:) = norm(v-u_true)*sqrt(h);
omega = [ 1; 1.2; 1.6; 1.95];

for j = 1:4
    w = omega(j);
for i = 1:iters
    v = WJac(A,f,v,w,1,relaxtype);
    ploterror(i+1,j) = norm(v-u_true)*sqrt(h);
end
end
figure(2);
set(gcf,'DefaultLineLineWidth',2);
plot(ploterror); legend(['omega = ',num2str(1)],['omega = ',num2str(1.2)],['omega = ',num2str(1.6)],['omega = ',num2str(1.95)]);
xlabel('iteration'); ylabel('error in l2 norm');