function [y]=tridiag(a,b,c,f)
         n=length(f);
         s=zeros(n,1);
         y=s;
         w=a(1);
         y(1)=f(1)/w;
         for ii=2:n
             s(ii-1)=c(ii-1)/w;
             w=a(ii)-b(ii)*s(ii-1);
             y(ii)=(y(ii)-b(ii)*y(ii-1))/w;
         end
         for ii=n-1:-1:1
             y(ii)=y(ii)-s(ii)*y(ii+1);
         end
end
