%% 1
clear all
close all
clc

L=pi;
N=1000;
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


e=abs((U-u)./U);
error=zeros(1,N);
for ii=1:N
    error(ii)=sum(e(2:N+1,ii+1))/N;
end
plot(t(2:N+1)',error)
%% 2


k=1;
w=0.001/(dt);

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
        F(jj,ii)=(w*cos(w*t(ii))+D*k^2*sin(w*t(ii)))*cos(k*x(jj));
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


e=abs((U-u)./U);
error=zeros(1,N);
for ii=1:N
    error(ii)=sum(e(2:N+1,ii+1))/N;
end
plot(t(2:N+1)',error)
