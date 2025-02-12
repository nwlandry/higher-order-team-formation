function log_likelihood = LL_2param(a,b,I0_collab,I0_no_collab)

% MLE Estimate
% Null model probability proportional to a*s_0+b

sumtest = ones(size(I0_collab,1),1);
for i=1: size(I0_collab,1)
    S = I0_collab(i);
    if S==0
        sumtest(i)=log(b);
    else
        sumtest(i) = log(a*S+b);
    end
end
mysum = nansum(sumtest);

sumtest2 = ones(size(I0_no_collab,1),1);
for i=1:size(I0_no_collab,1)
    S = I0_no_collab(i);
    if S==0
        sumtest2(i)=log(1-b);
    else
        sumtest2(i) = log(1-a*S-b);
    end
end
mysum2 = nansum(sumtest2);

log_likelihood = -(mysum + mysum2);