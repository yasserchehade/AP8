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
x=zeros(N,1);
for ii=1:N
    x(ii)=ii*dx;
end

t=zeros(N,1);
for ii=1:N
    t(ii)=ii*dt;
end

u=zeros(N,N+1);
u_0_t=0;
u_L_t=0;
u(:,1)=sin(k*x);
U=zeros(N);
for ii=1:N
    for jj=1:N
        U(ii,jj)=exp(-D*k^2*jj*dt)*sin(k*ii*dx);
    end
end

r=D*dt/dx^2;
A_n=2*(1-r)*eye(N)+r*diag(ones(N-1,1),1)+r*diag(ones(N-1,1),-1);
A_n_1=2*(1+r)*eye(N)-r*diag(ones(N-1,1),1)-r*diag(ones(N-1,1),-1);
for jj=2:N+1
    u(:,jj)=A_n_1\(A_n*u(:,jj-1));
end

u=u(:,2:N+1);

e=abs(U-u)./abs(U);
error=zeros(1,N);
for ii=1:N
    error(ii)=sum(e(:,ii))/N;
end
plot(t',error)
