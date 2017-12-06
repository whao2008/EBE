clear all

%Example 2 in the manuscript 
%True lambda
lambda = [2 16 24 96 -256 -1024];
global x w cx Tol
%Set the grid points. Here we use the uniform grid points
%Sparse grid can be also used, please see the following paper for more
%details.
%       Likelihood Approximation by Numerical Integration on Sparse Grids 
%       by Florian Heiss & Viktor Winschel
%       http://www.sparse-grids.de/
N=500;
x=(-1:2/N:1)';
w=ones(size(x))/length(x);
%Set the Tolerance for Newton's method
Tol=1e-7;
%Generate the PDF
n = length(lambda);
temp = lambda(1)*x;
for i=2:n
    temp = temp + lambda(i)*x.^i;
end
temp = exp(temp);
Z = w'*(temp);
p = temp/Z;

%Compute the moments
for i=1:n
    f(i) = w'*(x.^i.*p);
end

index=[1:6];
cx=generate_basis_matrix(x,index);

tic
%EBE
lambda_est=EMP_Newton_EBE_md(index,f);
%fsolve
%lambda_est=fsolve(@(lambda)
%nonlinear_fun_md_matrix(lambda,index,f),zeros(size(f))',...
%optimoptions(@fsolve,'SpecifyObjectiveGradient',true));
%Newton
%lambda_est=Newton_method(@(lambda) nonlinear_fun_md_matrix(lambda,index,f),zeros(size(f))');
toc


