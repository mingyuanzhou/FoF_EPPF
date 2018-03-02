function m= nk_to_m(n_k)
n_k=n_k(:)';
m=zeros(1,max(n_k));
for i=unique(n_k)
    m(i)=nnz(n_k==i);
end
  