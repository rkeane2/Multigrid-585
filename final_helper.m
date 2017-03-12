clear all;
figure(1)
n = 64;
h = 1/n;
probtype = 0;
x = (1:(n-1))*h; x = x';x = [0;x;1];
subplot(2,1,1);
u = Vcycle_backup(n,4,0*ones(n-1,1),probtype,3,3,2/3,0,1,10);
u = [1/2;u;1];
plot(x,u);

subplot(2,1,2); 
u_exact = exact(n,probtype);
plot(x,u_exact);

error = norm(u-u_exact)*sqrt(h)
error = norm(u-u_exact,inf);