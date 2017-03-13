function I = interpop(n,n2,restric,probtype,u)

% if probtype >2 %2-d case %NOTE: fix
%     n = n^2-2*n+2;
%     n2 = (n2^2)-2*n2+2;
%     
% end
if probtype <3
if restric ==0 %injection from n to n2
    I = zeros(n2-1,n-1); 
    for i=1:n2-1
        I(i,2*i) = 1;
    end
    I = I*u;
    return
end

I = zeros(n-1,n2-1); 
for i=1:n2-1
    I(2*i-1:2*i+1,i) = [1/2;1;1/2]; 
end
 
if restric ==1
    I = (1/2)*I'; %restriction full weighting
    I = I*u;
    return
end

if restric ==2
    
    
%     if probtype ==0
%         l_b = 1/2;
% r_b = 1;
%     end
% 
% if probtype ==1
%         l_b = 1/2;
% r_b = -1/2;
% end
% if probtype ==2
%         l_b = 1;
% r_b = 3;
% end
% 
%     I = [zeros(n-1,1),I,zeros(n-1,1)];
%     I(1,1) = 1/2; I(end,end) = 1/2;
%     u = [l_b;u;r_b];
    I = I*u;
end
    


end
%interp operator on page 40 

%%
if probtype>2 %2d interpolation/restriction
    if restric ==0 %injection 
            n = n^2-2*n+2;
    n2 = (n2^2)-2*n2+2;
    I = zeros(n2-1,n-1);
    for i=1:n2-1
        I(i,2*i) = 1;
    end
    I = I*u;
    return
    end
    
    if restric ==1 %full weighting restriction 
        n = sqrt(length(u))+1;
        u = reshape(u,n-1,n-1);
        I = zeros(n/2-1,n/2-1);
        for i=1:n/2-1
            for j=1:n/2-1
                I(i,j) = 1/16*(u(2*i-1,2*j-1)+u(2*i-1,2*j+1)+u(2*i+1,2*j+1)+...
                    u(2*i+1,2*j-1)+2*(u(2*i,2*j-1)+u(2*i,2*j+1)+u(2*i-1,2*j)+...
                    u(2*i+1,2*j))+4*u(2*i,2*j));
            end
        end
        I = I(:);
    end
    
    if restric ==2 %interpolation
        n = sqrt(length(u))+1;
        u = reshape(u,n-1,n-1);
        I = zeros(2*n+1,2*n+1);
        %if probtype >2
            u = [zeros(1,n+1);zeros(n-1,1),u,zeros(n-1,1);zeros(1,n+1)];
        %end
        for i=1:n
            for j=1:n
                I(2*i-1,2*j-1) = u(i,j);
                I(2*i,2*j-1) = 1/2*(u(i,j)+u(i+1,j));
                I(2*i-1,2*j) = 1/2*(u(i,j)+u(i,j+1));
                I(2*i,2*j) = 1/4*(u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1));
            end
        end
        I = I(2:end-1,2:end-1); 
        %u = u(2:end-1.2:end-1); u = u(:);
        %I = I*u;
%         I(1,1) = 1/4*u(1,1);
%         I(1,end) = 1/4*u(1,end); 
%         I(end,1) = 1/4*u(end,1); 
%         I(end,end) = 1/4*u(end,end);       
        %         for j=1:n-2
%             I(1,2*j) = 1/2*(u(1,j));
%             I(1,2*j+1) = 1/4*(u(1,j)+u(1,j+1));
%             I(2*j,1) = 1/2*(u(j,1));
%             I(2*j+1,1) = 1/4*(u(j,1)+u(j+1,1));
%             I(end,2*j) = 1/2*(u(end,j));
%             I(end,2*j+1) = 1/4*(u(end,j)+u(end,j+1));
%             I(2*j,end) = 1/2*(u(j,end));
%             I(2*j+1,end) = 1/4*(u(j,end)+u(j+1,end));
%         end
        I = I(:);
    end
end
        
        
                
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
