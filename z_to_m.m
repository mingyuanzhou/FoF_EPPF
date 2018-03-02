function m= z_to_m(z)
n_k=sparse(1,z,1);
m=zeros(1,max(n_k));
for i=1:max(n_k)
    m(i)=nnz(n_k==i);
end
  