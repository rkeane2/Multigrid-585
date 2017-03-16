clear all;

%% test solution 1-d
figure(1); clf;
n =256;
%v = linspace(1,3,n-1)';
 probtype = 1;
 levels =3;
 w = 1.2; 
 relaxtype = 1;
 resttype = 1;
 iters = 5;
 v = ones(n-1,1); 
  %v = exact(n,probtype);
  %v = v(2:end-1);

% w = 1;
% relaxtype = 1;
 
% x = (1:(n-1))*h; x = x';x = [0;x;1];
subplot(2,1,1);
 %ylim([1 6])

%subplot(3,1,2); 
[u_exact,x] = exact(n,probtype);
hold on; plot(x,u_exact); 
[u] = Vcycle(n,levels,v,probtype,3,3,w,relaxtype,resttype,iters,1);
% u = [1/2;u;1];
plot(x,u,'r-'); title('computed and true solutions');

h = 1/n;
error_l2 = norm(u-u_exact)*sqrt(h)
residual = norm(matrix(n,probtype)*u(2:end-1)-rhs(n,probtype))*sqrt(h);
%error_inf = norm(u-u_exact,inf)

subplot(2,1,2);

error = u-u_exact;
plot(x,error); title(['||error|| is ',num2str(error_l2),' after ',num2str(iters),' iterations'])

%% test sol 2-d 
figure(1)
n =64;

%v = linspace(1,3,n-1)';
 probtype = 5;
 levels = 2;
 w = 1.3; 
 relaxtype = 1;
 resttype = 1;
 iters = 1;
 v = 1*ones((n-1)^2,1); 

subplot(3,1,2); 

[u_exact,x] = exact(n,probtype); 

surf(x(1:n+1,1:n+1),x(1:n+1,n+2:2*n+2),u_exact); shading interp;title('true solution'); view(-90,90);

 subplot(3,1,1);

[u] = Vcycle(n,levels,v,probtype,3,3,w,relaxtype,resttype,iters,1);
surf(x(1:n+1,1:n+1),x(1:n+1,n+2:2*n+2),u); shading interp; title('computed solution');view(-90,90);

h = 1/n;
error_l2 = norm(u(:)-u_exact(:))*h
error_inf = norm(u(:)-u_exact(:),inf)
subplot(3,1,3); 
error = u-u_exact;
surf(x(1:n+1,1:n+1),x(1:n+1,n+2:2*n+2),error); shading interp; title(['||error|| is ',num2str(error_l2),' after ',num2str(iters),' iterations'])
view(-95,38);
%% test nonlinear problem

figure(1)
n =64; %can show how this changes the error
 probtype = 6;
 levels = 2;
 w = 1; 
 relaxtype = 0;
 resttype = 1;
 iters = 5;
 %v = ones((n-1)^2,1); %won't converge with this guess
v = zeros((n-1)^2,1);

subplot(3,1,2);

[u_exact,x] = exact(n,probtype); 

surf(x(1:n+1,1:n+1),x(1:n+1,n+2:2*n+2),u_exact); shading interp;view(-90,90);

 subplot(3,1,1);
 u = NMG(n,levels,v,3,3,w,relaxtype,resttype,iters);
 surf(x(1:n+1,1:n+1),x(1:n+1,n+2:2*n+2),u); shading interp;view(-90,90);
 
 h = 1/n;
error_l2 = norm(u(:)-u_exact(:))*h
error_inf = norm(u(:)-u_exact(:),inf)
subplot(3,1,3); 
error = u-u_exact;
surf(x(1:n+1,1:n+1),x(1:n+1,n+2:2*n+2),error); shading interp;title(['||error|| is ',num2str(error_l2),' after ',num2str(iters),' iterations'])

view(-95,38);


















%%

A = matrix(256,2); 
f = rhs(256,2); 
v = zeros(255,1);
w=.5;

test1 = WJac(A,f,v,w,10,0);
test2 = norm(A*test1-f)