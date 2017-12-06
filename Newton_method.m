function [x,flag]=Newton_method(fun_handle,x)
global Tol;
%Tol=1e-7;
Display=1;
MaxIter=100;
ct=1;
flag=0;
while 1
    [f,df]=feval(fun_handle,x);
    if any(any(isnan(df) | isinf(df)))
        flag=-1;
        return;
    end
    if norm(f)<Tol
        flag=1;
        break
    end
    dx=-df\f;
    if norm(dx)<Tol
        flag=2;
        break;
    end
    x=x+dx;
    ct=ct+1;
    if ct>MaxIter
        flag=-2;
        break;
    end
end
global gct
gct=gct+ct;
if Display
    if flag==1
    %    disp(['    Equation solved after ' num2str(ct) ' iterations!']);
    end
    if flag==2
        disp(['    dx is too small after ' num2str(ct) ...
            ' iterations, residual is ' num2str(norm(f))]);
    end
    if flag==-2
        disp(['    failed, reach maximium iterations...']);
    end    
end
