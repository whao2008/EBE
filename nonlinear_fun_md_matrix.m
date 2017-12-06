function [ff,df,dp]=nonlinear_fun_md_matrix(lambda,index,f,lambda2)
global cx x w
%load spgpt
d=size(index,1);
n=length(f);

if nargin>3
    lambda=[lambda;lambda2];
    m=n+1;
else
    m=n;
end

va=exp(cx(:,1:m)*lambda(1:m));

if nargout>1
    df=zeros(n,n);
end
if nargout>2
    dp=zeros(n,1);
end
ff=zeros(n,1);
for i=1:n
    ff(i)=(w'*(cx(:,i).*va))/(w'*va)-f(i);
    if nargout>1
        for j=1:n
            df(i,j)=((w'*(cx(:,i).*cx(:,j).*va))*(w'*va)-...
                (w'*(cx(:,i).*va))*(w'*(cx(:,j).*va)))/(w'*va)^2;
        end
    end
    if nargout>2
        temp1=1;
        for j1 = 1:d
            temp1 = cx(:,j).*x(:,j1).^index(j1,m);
        end
        dp(i)=((w'*(cx(:,i).*cx(:,j).*va))*(w'*va)-...
            (w'*(cx(:,i).*va))*(w'*(cx(:,j).*va)))/(w'*va)^2;
    end
end
% f=f(end:-1:1);
% if nargout>1
%     df=df(end:-1:1,:);
%     if nargout>2
%         dp=dp(end:-1:1);
%     end
% end


