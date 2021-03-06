
%note: the relevant quantities are residualhold, residualratio, errorhold,
%and errorratio. We denote the residuals for V-cycle as above, and add a
%"f" at the end for FMG, and add a "n" for NMG.


%% 1-D case V cycle

n =256;

 probtype = 2;
 levels =3;
 w = 1.1; 
 relaxtype = 1;
 resttype = 1;
 iters = 15;
 v = ones(n-1,1);
 A = matrix(n,probtype); 
 f = rhs(n,probtype);
 
 [u_exact,x] = exact(n,probtype);
 residualhold = zeros(iters,1); 
 errorhold = zeros(iters,1); 
 residualratio = zeros(iters-1,1); 
 errorratio = zeros(iters-1,1);
 
 for i = 1:iters
 [u] = Vcycle(n,levels,v,probtype,3,3,w,relaxtype,resttype,1,1);
 v = u(2:end-1);
 residual = norm(A*v-f)*sqrt(1/n);
 error = norm(u_exact-u)*sqrt(1/n) ;
 residualhold(i) = residual;
 errorhold(i) = error;
 if i==1
     continue
 end
 errorratio(i-1,1) = errorhold(i)/errorhold(i-1);
 residualratio(i-1,1) = residualhold(i)/residualhold(i-1);
 end
 residualratio;
 
 
 %% 2-d case V-cycle
 
 n =64;

 %n =128;

%v = linspace(1,3,n-1)';
 probtype = 4;
 levels = 3;
 w = 1.3; 
 relaxtype = 1;
 resttype = 1;
 iters = 12;
 v = 1*ones((n-1)^2,1); 
 %v = ones((n-1)^2,1);
 A = matrix(n,probtype); 
 f = rhs(n,probtype);
 
 [u_exact,x] = exact(n,probtype);
 residualhold = zeros(iters,1); 
 errorhold = zeros(iters,1); 
 residualratio = zeros(iters-1,1); 
 errorratio = zeros(iters-1,1);
 
 for i = 1:iters
 [u] = Vcycle(n,levels,v,probtype,3,3,w,relaxtype,resttype,1,1);
 v = u(2:end-1,2:end-1);
 v = v(:);
 residual = norm(A*v-f)*(1/n);
 error = norm(u_exact-u)*(1/n) ;
 residualhold(i) = residual;
 errorhold(i) = error;
 if i==1
     continue
 end
 errorratio(i-1,1) = errorhold(i)/errorhold(i-1);
 residualratio(i-1,1) = residualhold(i)/residualhold(i-1);
 end
 residualratio;
 
 %% 1-D case FMG 
 
 n =256;

 probtype = 2;
 levels =3;
 w = 1.1; 
 relaxtype = 1;
 resttype = 1;
 iters = 15;
 v = ones(n-1,1);
 A = matrix(n,probtype); 
 f = rhs(n,probtype);
 
 [u_exact,x] = exact(n,probtype);
 residualholdf = zeros(iters,1); 
 errorholdf = zeros(iters,1); 
 residualratiof = zeros(iters-1,1); 
 errorratiof = zeros(iters-1,1);
 
 for i = 1:iters
 [u] = FMG(n,levels,probtype,3,3,w,relaxtype,resttype,i,1);
 v = u(2:end-1);
 residual = norm(A*v-f)*sqrt(1/n);
 error = norm(u_exact-u)*sqrt(1/n) ;
 residualholdf(i) = residual;
 errorholdf(i) = error;
 if i==1
     continue
 end
 errorratiof(i-1,1) = errorholdf(i)/errorholdf(i-1);
 residualratiof(i-1,1) = residualholdf(i)/residualholdf(i-1);
 end
 residualratio;
 
 %% 2-D FMG 
 
  n =64;

 probtype = 4;
 levels =3;
 w = 1.3; 
 relaxtype = 1;
 resttype = 1;
 iters = 10;
 v = ones((n-1)^2,1);
 A = matrix(n,probtype); 
 f = rhs(n,probtype);
 
 [u_exact,x] = exact(n,probtype);
 residualholdf = zeros(iters,1); 
 errorholdf = zeros(iters,1); 
 residualratiof = zeros(iters-1,1); 
 errorratiof = zeros(iters-1,1);
 
 for i = 1:iters
 [u] = FMG(n,levels,probtype,3,3,w,relaxtype,resttype,i,1);
 v = u(2:end-1,2:end-1);
 v = v(:);
 residual = norm(A*v-f)*(1/n);
 error = norm(u_exact-u)*(1/n) ;
 residualholdf(i) = residual;
 errorholdf(i) = error;
 if i==1
     continue
 end
 errorratiof(i-1,1) = errorholdf(i)/errorholdf(i-1);
 residualratiof(i-1,1) = residualholdf(i)/residualholdf(i-1);
 end
 residualratio;
 
 
 %% 2-D nonlinear
 
 n =64;

 probtype = 6;
 levels =2;
 w = .5; 
 relaxtype = 0;
 resttype = 1;
 iters = 10;
 v = zeros((n-1)^2,1);
 
 [u_exact,x] = exact(n,probtype);
 residualholdn = zeros(iters,1); 
 errorholdn = zeros(iters,1); 
 residualration = zeros(iters-1,1); 
 errorration = zeros(iters-1,1);
  for i = 1:iters
 [u] = NMG(n,levels,v,3,3,w,relaxtype,resttype,1);
 v = u(2:end-1,2:end-1);
 v = v(:);
 residual = norm(rhs(n,probtype,v))*(1/n);
 error = norm(u_exact-u)*(1/n) ;
 residualholdn(i) = residual;
 errorholdn(i) = error;
 if i==1
     continue
 end
 errorration(i-1,1) = errorholdn(i)/errorholdn(i-1);
 residualration(i-1,1) = residualholdn(i)/residualholdn(i-1);
 end
 