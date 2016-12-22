clear
close all

format long
a=-.02;
xi = 0;
xf = 1;

ti = 0;
tf = .05;


%calculating exact 
N = 200;
    
    h = (xf-xi)/N;
    delt=.00005;
    lambda=delt/(h^2);
    M = fix((tf-ti)/delt);
    u0 = zeros(N-1,1);
    X=zeros(N-1,1);
    Y = zeros(N-1,1);
    u1 = zeros(N-1,1);  
    F = zeros(N-1,1);
    uf1 = zeros(N-1,N-1);
    
    for i=1:N-1
        X(i)=xi+i*h;
        u0(i) = sin(pi*X(i));
    end

    t=0;
    count1 = 0;
    
    for j = 1:M
    
        tol = .0000001;
        err = 1;
        
        
    
        while(tol<err)
            F(1) = u1(1)*(1 + 2*lambda + a*delt) - (a+1)*delt*(u1(1))^2 + delt*(u1(1))^3  - lambda*u1(2) - u0(1);
            F(N-1) = u1(N-1)*(1 + 2*lambda + a*delt) - (a+1)*delt*(u1(N-1))^2 + delt*(u1(N-1))^3 - lambda*u1(N-2) - u0(N-1);
   
        for i = 2:N-2
            F(i) = u1(i)*(1 + 2*lambda + a*delt) - (a+1)*delt*(u1(i))^2 + delt*(u1(i))^3 - lambda*u1(i-1) - lambda*u1(i+1) - u0(i);
        end
        
        J = zeros(N-1,N-1);
        
        for i = 1:N-2
            J(i,i) = (1 + 2*lambda + a*delt) - 2*(a+1)*delt*u1(i) + 3*delt*(u1(i)^2);
            J(i+1,i) = -lambda;
            J(i,i+1) = -lambda;
        end
        
        J(N-1,N-1) = (1 + 2*lambda + a*delt) - 2*(a+1)*delt*u1(N-1) + 3*delt*(u1(N-1)^2);
            
        p = J\(-F);        
        
        u1J = u1 + p;
        err = max(abs(p));
        u1 = u1J;
         
        count1 = count1 +1;
        end
        for o=1:N-1
            uf1(j,o) = u1(o);
        end
        
        
    u0 = u1;
    t = t + delt;
    end
    exact = u1;
    
   
% calculating sol for more refined grid  
N = 10;
exact1 = zeros(3,39);
for k=1:3
    h = (xf-xi)/N;
    L(k)=h;
    delt=0.00005;
    lambda=delt/(h^2);
    M = fix((tf-ti)/delt);
    u0 = zeros(N-1,1);
    X=zeros(N-1,1);
    Y = zeros(N-1,1);
    u1 = zeros(N-1,1);  
    F = zeros(N-1,1);
    uf = zeros(N-1,N-1);
    
    for i=1:N-1
        X(i)=xi+i*h;
        u0(i) = sin(pi*X(i));
    end

    t=0;
    count = 0;
    for i=1:M
        Y(i) = ti + i*delt;
    end
    
    for j = 1:M
    
        tol = .0000001;
        err = 1;
        
        
    
        while(tol<err)
            F(1) = u1(1)*(1 + 2*lambda + a*delt) - (a+1)*delt*(u1(1))^2 + delt*(u1(1))^3  - lambda*u1(2) - u0(1);
            F(N-1) = u1(N-1)*(1 + 2*lambda + a*delt) - (a+1)*delt*(u1(N-1))^2 + delt*(u1(N-1))^3 - lambda*u1(N-2) - u0(N-1);
        
            for i = 2:N-2
                F(i) = u1(i)*(1 + 2*lambda + a*delt) - (a+1)*delt*(u1(i))^2 + delt*(u1(i))^3 - lambda*u1(i-1) - lambda*u1(i+1) - u0(i);
            end
        
            J = zeros(N-1,N-1);
        
            for i = 1:N-2
                J(i,i) = (1 + 2*lambda + a*delt) - 2*(a+1)*delt*u1(i) + 3*delt*(u1(i)^2);
                J(i+1,i) = -lambda;
                J(i,i+1) = -lambda;
            end
        
            J(N-1,N-1) = (1 + 2*lambda + a*delt) - 2*(a+1)*delt*u1(N-1) + 3*delt*(u1(N-1)^2);
            
            p = J\(-F);
        
         
            u1J = u1 + p;
            err = max(abs(p));
            u1 = u1J;
         
            count = count +1;
        end
    
        for o=1:N-1
            uf(j,o) = u1(o);
        end
    u0 = u1;
    t = t + delt;
    end   
    for i=1:N-1
         exact1(k,i) = exact((i*(200/N)));
    end
    
    error(k) = max(abs(exact1(k,1:N-1)'- u1));
    
    N=2*N;
end
    
  for i= 1:2
       order(i)= log(error(i)/error(i+1))/log(L(i)/L(i+1))
  end
  
    uf2 = zeros(1000,39);
    for j=1:M
        for o=1:39
           uf2(j,o) = uf1(j,o*(200/40));
        end
    end
figure(1)
surf(X,Y,uf);
shading interp
colorbar
xlabel('X axis (Space)')
ylabel('Y axis (Time)')
zlabel('U(x,t)')
title('approximate')
figure(2)
surf(X,Y,uf2);
shading interp
colorbar
xlabel('X axis (Space)')
ylabel('Y axis (Time)')
zlabel('U(x,t)')
title('exact')

ufy = uf(:,20);
uf2y = uf2(:,20);
ufy1 = uf(:,4);
uf2y1 = uf2(:,4);

figure(3)
plot(Y,ufy,Y,ufy1);
xlabel('Time')
ylabel('U(.5,t)')
title('approximate (red:a=.1; blue:a=.5)')

figure(4)
plot(Y,uf2y,'.',Y,uf2y1,'.');
xlabel('Time')
ylabel('U(a,t)')
title('exact (red:a=.1; blue:a=.5)')

figure(5)
plot(X,u1,X,exact1(3,:),'o');
title('approximate:"-"  exact:"o"')
