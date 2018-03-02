function [z,n_k] = m_i_to_z(m_j,jj)
%Convert {m_1,m_2,...} to (z_1,z_2,...,z_n) and (n_1,n_2,...,n_l)
Lenx = sum(m_j.*jj);
z = zeros(1,Lenx);
count=0;
label=0;
for ii=1:length(jj)
    j=jj(ii);
 
    for i=1:m_j(ii)
        label=label+1;
        z(count+(1:j)) = label;
        count=count+j;
    end
end

n_k=sparse(z,1,1);
% 
% n_k=zeros(length(unique(z)),1);
% count=0;
% for k=unique(z)
%     count=count+1;
%     n_k(count)=nnz(z==k);
% end