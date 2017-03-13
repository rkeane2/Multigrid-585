function u = NMG(n,levels,v,mu1,mu2,w,relaxtype,resttype,iters)

for i = 1:iters
     %v = v+JM(v)\rhs(n,4,v);
    cor = Vcycle(n,levels,v,5,mu1,mu2,w,relaxtype,resttype,1,0);
    v = v+cor;
end
u = reshape(v,n-1,n-1);
u = [zeros(1,n+1);zeros(n-1,1),u,zeros(n-1,1);zeros(1,n+1)];