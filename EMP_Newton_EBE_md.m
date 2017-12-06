function [lambda_success,index,moments]=EMP_Newton_EBE_md(index,moments)
n=length(moments);
Tol0=5e-1;
lambda=0;
global cx
i=1;
while 1
    Tol=Tol0;
    while 1
        lambda0=lambda;
        [lambda,flag]=Newton_method(@(lambda) nonlinear_fun_md_matrix(lambda,...
            index,moments(1:i)),lambda);
        if flag>=0
            lambda_success=lambda;
            break
        else            
            lambda=lambda0;
            Tol=Tol/10;
            [lambda,flag]=EMP_Newton_EBE_homotopy(index,moments,lambda(end),lambda(1:end-1),Tol,i-1);
            if flag
                if flag==2
                    lambda=lambda0;
                end
                disp(['remove ' num2str(i) '-th moment'])
                if i<size(index,2)
                    index(:,i)=[];
                    moments(i)=[];
                    cx(:,i)=[];
                    lambda=lambda_success;
                    i=i-1;
                    break
                end
                if i>size(index,2)
                    return
                end
            end
        end
    end
    if i==size(index,2)
        return
    end
    lambda1{i}=lambda;
    %save solsKS4d lambda1
    %disp([num2str(i) '-th: ' num2str(lambda)])
    Tol=Tol0;
    while 1
        lambda0=lambda;
        [lambda,flag]=EMP_Newton_EBE_homotopy(index,moments,0,lambda0,Tol,i);
        if flag
            if i<size(index,2)
                disp(['remove ' num2str(i+1) '-th moment'])
                index(:,i+1)=[];
                moments(i+1)=[];
                cx(:,i+1)=[];
                lambda=lambda0;
            end
            if i==size(index,2)
                return
            end
        else
            break
        end
    end
    if i==56
    end
    i=i+1
end
