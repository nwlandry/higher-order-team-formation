function log_likelihood = LL1(p,collab_values,no_collab_values)

% MLE Estimate
%For the null model where all of the participants have a constant
%probability of collaborating. 

sumtest = ones(size(collab_values,1),1);
for i=1: size(collab_values,1)
    sumtest(i) = log(p);
end
mysum = sum(sumtest);

sumtest2 = ones(size(no_collab_values,1),1);
for i=1:size(no_collab_values,1)
    sumtest2(i) = log(1-p);
end
mysum2 = sum(sumtest2);

log_likelihood = -(mysum + mysum2);