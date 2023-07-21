function [h,p]=realdata_bootstrap(d1,d2,n,alpha)
%d1,d2: input data (PSTHs)
%n: number of bootstrap instances
%alpha: significance level

real_diff=nanmean(d1)-nanmean(d2);

values=[d1;d2];
ids=[ones(size(d1));2*ones(size(d2))];

b_vec=[];
parfor i=1:n
    ids_b=ids(randperm(length(ids)));
    b_vec(i)=nanmean(values(find(ids_b==1)))-nanmean(values(find(ids_b==2)));
end

p=1-length(find(b_vec<real_diff))/length(b_vec);

if p<alpha
    h=1;
else
    h=0;
end
