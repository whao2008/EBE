function [lambda,flag_hom]=EMP_Newton_EBE_homotopy(index,moments,lambda2,lambda,Tol,i)
ii=0;
MaxIter=10;
flag_hom=0;
while 1
    [fi,df,dp]=nonlinear_fun_md_matrix([lambda;lambda2],index,moments(1:i+1));
    fi=fi(end);
    if abs(fi)<Tol
        lambda=[lambda;lambda2];
        break
    end
    ii=ii+1;
    if ii>MaxIter
        lambda=[lambda;lambda2];
        break
    end
    disp([num2str(ii) ':' num2str(lambda2) ',' num2str(fi)]);
    [dfi,dp1]=dfnonlinear_fun_md_matrix([lambda;lambda2],index,moments(1:i+1),i+1);
    dlambda2_final=-((-df(1:i,1:i)\dp(1:i))'*dp1+dfi)\fi;
    dlambda2=[0;dlambda2_final];
    dlambda=-df(1:i,1:i)\dp(1:i);
    while 1
        [lambda0,flag]=Newton_method(@(lambda) nonlinear_fun_md_matrix(lambda,...
            index,moments(1:i),lambda2+dlambda2(end)),lambda+dlambda*dlambda2(end));
        %         for j=1:i
        %             subplot(1,i,j)
        %             hold on
        %             plot(lambda2+dlambda2(end),lambda0(j),'r*')
        %         end
        disp(['   ' num2str(flag) ':' num2str(dlambda2(end))])
        if flag==1
            lambda=lambda0;
            lambda2=lambda2+dlambda2(end);
            if abs(dlambda2(end)-dlambda2_final)<1e-13
                break
            else
                if abs(dlambda2(end))<0.01*abs(dlambda2_final)
                    lambda=[lambda;lambda2];
                    flag_hom=1;
                    return
                end
                dlambda2_final=dlambda2_final-dlambda2(end);
                [fi,df,dp]=nonlinear_fun_md_matrix([lambda;lambda2],index,moments(1:i+1));
                %                 hold on
                %                 plot(lambda2,fi(end),'r*')
                if abs(fi(end))<Tol
                    break
                end
                dlambda=-df(1:i,1:i)\dp(1:i);
                dlambda2=[0;sign(dlambda2(end))*min(abs(dlambda2(end)),abs(dlambda2_final))];
            end
        else
            dlambda2=[dlambda2(1:end-1);(dlambda2(end)+dlambda2(end-1))/2];
            if abs(dlambda2(end))<0.1*abs(dlambda2_final)
                flag_hom=2;
                return
            end
        end
    end
    %     if cur_Iter>MaxIter
    %         disp('max iter is reached')
    %         lambda=[lambda;lambda2];
    %         break
    %     end
end







