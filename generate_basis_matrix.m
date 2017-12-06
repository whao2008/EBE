function cx=generate_basis_matrix(x,index)
d=size(x,2);
n=size(x,1);
cx=ones(n,size(index,2));
for i=1:size(index,2)
    for j = 1:d
        cx(:,i) = cx(:,i).*x(:,j).^index(j,i);
    end
end