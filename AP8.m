%% 1
clear all
close all
clc

L=pi;
N=100;
T=10;
dx=L/(N+1);
dt=T/(N+1);
k=1;
D=0.1;
F=zeros(N+2);
x=zeros(N+2,1);
for ii=1:N+1
    x(ii+1)=ii*dx;
end

t=zeros(1,N+2);
for ii=1:N+1
    t(ii+1)=ii*dt;
end

u=zeros(N+2);
u(:,1)=sin(k*x);
U=zeros(N);
for ii=1:N+2
    for jj=1:N+2
        U(ii,jj)=exp(-D*k^2*t(jj))*sin(k*x(ii));
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
%% graphing 1
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
u=zeros(N+2);
u(1,:)=sin(w*t);
u(N+2,:)=sin(w*t)*cos(k*L);
U=zeros(N);
for ii=1:N+2
    for jj=1:N+2
        U(ii,jj)=sin(w*t(jj))*cos(k*x(ii));
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
%% graphing 2
for a=[1/4 1/2 3/4]
    j=round((N+1)*a);
    figure
    plot(x,u(:,j),x,U(:,j))
    legend('calculated','exact')
    title(['plot comaprison for t=' num2str(T*a) ' of calculated and exact' ])
    xlabel('x')
    ylabel('u(x,t)')
  
end


%%
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


    e=abs((U-u)./U);
    error(z)=sum(e(2:N+1,N+2))/N;
end
plot(w,error)
title(['error vs w'])
xlabel('w')
ylabel('error')
