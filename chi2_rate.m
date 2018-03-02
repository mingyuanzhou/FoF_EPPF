function xi = chi2_rate(M, Rate,ii)

estimated = Rate;
if nargin<3
    for ii=length(Rate):-1:1
        if estimated(ii)>=1
            break
        end
    end
    ii=ii+1;
end

%ii=5;

%ii = find(estimated<4);
estimated(ii)=sum(Rate(ii:end));
estimated = estimated(1:ii);


xi = [M(1:ii-1),sum(M(ii:end))];
xi = sum((xi-estimated).^2./estimated);