function [ff,df]=dfnonlinear_fun_md_matrix(lambda,index,f,ii)
global cx x w 
%load spgpt
m=length(f);
va=exp(cx(:,1:m)*lambda(1:m));
ff=((w'*(cx(:,ii).^2.*va))*(w'*va)-(w'*(cx(:,ii).*va))*(w'*(cx(:,ii).*va)))/(w'*va)^2;

if nargout>1
    df=zeros(ii-1,1);
    for i=1:ii-1
        df(i)=((w'*(cx(:,ii).*cx(:,i).*va))*(w'*va)-...
            (w'*(cx(:,ii).*va))*(w'*(cx(:,i).*va)))/(w'*va)^2;
    end
end
end

