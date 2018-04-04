%% 1
clear all
close all
clc
%defining parameters
L=pi;
N=100;
T=10;
dx=L/(N+1);
dt=T/(N+1);
k=1;
D=0.1;
%creating a matrix to fill up the known vectors
F=zeros(N+2);
%defining the x axis
x=zeros(N+2,1);
for ii=1:N+1
    x(ii+1)=ii*dx;
end
%defining the t axis
t=zeros(1,N+2);
for ii=1:N+1
    t(ii+1)=ii*dt;
end
%defining a matrix for the solution
u=zeros(N+2);
%column 1 and N+2 are at boundaris 0 and T respectively 
%row 1 and N+2 are at boundaries 0 and L respectively 
%inputing boundary condition
u(:,1)=sin(k*x);

%defining the exact solution
U=zeros(N);
for ii=1:N+2
    for jj=1:N+2
        U(ii,jj)=exp(-D*k^2*t(jj))*sin(k*x(ii));
    end
end

%defining our lambda
r=D*dt/dx^2;
%defining our RHS tridiagonal matrix
A_n=2*(1-r)*eye(N)+r*diag(ones(N-1,1),1)+r*diag(ones(N-1,1),-1);
%defining the parameters of the LHS tridiagonal matrix
a=2*(1+r)*ones(N,1);
b=-r*ones(N,1);
c=b;
%the loop is to do the calculation for each column respectively going from
%t to t+dt each time
for jj=2:N+2
    %solving the RHS
    f=A_n*u(2:N+1,jj-1)+2*dt*F(2:N+1,jj-1);
    %adding the boundary values 
    f(1,1)=f(1,1)+r*(u(1,jj)+u(1,jj-1));
    f(N,1)=f(N,1)+r*(u(N+2,jj)+u(N+2,jj-1));
    %solving the equation using the tridiagonal algorithm 
    u(2:N+1,jj)=tridiag(a,b,c,f);
end
%% graphing 1
%ploting the graphs 
for a=[1/4 1/2 3/4]
    j=round((N+1)*a);
    figure
    plot(x,u(:,j),x,U(:,j))
    legend('calculated','exact')
    title(['plot comaprison for t=' num2str(T*a) ' of calculated and exact' ])
    xlabel('x')
    ylabel('u(x,t)')
  
end


%% 2
%similarly defining parameters
k=1;
w=0.1/(dt);

x=zeros(N+2,1);
for ii=1:N+1
    x(ii+1)=ii*dx;
end

t=zeros(1,N+2);
for ii=1:N+1
    t(ii+1)=ii*dt;
end
F=zeros(N+2);
for ii=1:N+2
    for jj=1:N+2
        F(ii,jj)=(w*cos(w*t(jj))+D*(k^2)*sin(w*t(jj)))*cos(k*x(ii));
    end
end
%adding boundary counditions
u=zeros(N+2);
u(1,:)=sin(w*t);
u(N+2,:)=sin(w*t)*cos(k*L);
%geting the exact solution
U=zeros(N);
for ii=1:N+2
    for jj=1:N+2
        U(ii,jj)=sin(w*t(jj))*cos(k*x(ii));
    end
end
%defining the known parmeters for the RHS
r=D*dt/dx^2;
A_n=2*(1-r)*eye(N)+r*diag(ones(N-1,1),1)+r*diag(ones(N-1,1),-1);
%defining the parameters for the LHS
a=2*(1+r)*ones(N,1);
b=-r*ones(N,1);
c=b;
%solution
for jj=2:N+2
    %solving the RHS
    f=A_n*u(2:N+1,jj-1)+2*dt*F(2:N+1,jj-1);
    %adding the boundary counditions
    f(1,1)=f(1,1)+r*(u(1,jj)+u(1,jj-1));
    f(N,1)=f(N,1)+r*(u(N+2,jj)+u(N+2,jj-1));
    %solving the LHS
    u(2:N+1,jj)=tridiag(a,b,c,f);
 
end
%% graphing 2
%ploting the results
for a=[1/4 1/2 3/4]
    j=round((N+1)*a);
    figure
    plot(x,u(:,j),x,U(:,j))
    legend('calculated','exact')
    title(['plot comaprison for t=' num2str(T*a) ' of calculated and exact' ])
    xlabel('x')
    ylabel('u(x,t)')
  
end


%% effect of w
% solving the previous problem with diferent increasing w

F=zeros(N+2);
u=zeros(N+2);
U=zeros(N);
v=60;
error=zeros(1,v);
w=zeros(1,v);
for z=1:v
    w(z)=0.025*(z)/dt;

    for ii=1:N+2
        for jj=1:N+2
            F(jj,ii)=(w(z)*cos(w(z)*t(ii))+D*k^2*sin(w(z)*t(ii)))*cos(k*x(jj));
        end
    end

    u(1,:)=sin(w(z)*t);
    u(N+2,:)=sin(w(z)*t)*cos(k*L);

    for ii=1:N+2
        for jj=1:N+2
            U(ii,jj)=sin(w(z)*t(jj))*cos(k*x(ii));
        end
    end

    r=D*dt/dx^2;
    A_n=2*(1-r)*eye(N)+r*diag(ones(N-1,1),1)+r*diag(ones(N-1,1),-1);

    a=2*(1+r)*ones(N,1);
    b=-r*ones(N,1);
    c=b;
    for jj=2:N+2
        f=A_n*u(2:N+1,jj-1)+2*dt*F(2:N+1,jj-1);
        f(1,1)=f(1,1)+r*(u(1,jj)+u(1,jj-1));
        f(N,1)=f(N,1)+r*(u(N+2,jj)+u(N+2,jj-1));
        u(2:N+1,jj)=tridiag(a,b,c,f);
    end

%calculating the error
    e=abs((U-u)./U);
    error(z)=sum(e(2:N+1,N+2))/N;
end
%ploting the error with w
plot(w,error)
title(['error vs w'])
xlabel('w')
ylabel('error')
