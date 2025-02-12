function log_likelihood = logisticNegLogLike(params, allData, idxSubset)
    % params: [beta0, beta1, beta2, ..., betaK] for the chosen subset
    % X: table or matrix of *all* covariates
    % y: 0/1 vector
    % idxSubset: which columns of X to actually use

    % First param is intercept
    b = params(1);

    % Next ones are for the columns in idxSubset
    beta = params(2:end);

    % Collaboration
    allData_collab = allData(allData.collaborated==1,:);
    allData_noCollab = allData(allData.collaborated==0,:);

    X_collab = [allData_collab.sync1, allData_collab.async1, allData_collab.async2, allData_collab.async3];
    X_noCollab = [allData_noCollab.sync1, allData_noCollab.async1, allData_noCollab.async2, allData_noCollab.async3];

    % Pull out only the columns we need
    Xs_collab = X_collab(:, idxSubset);
    Xs_noCollab = X_noCollab(:, idxSubset);

    % Compute linear predictor = beta0 + sum_j( beta_j * Xs_j )
    linPred_collab = log(b + Xs_collab*beta');
    linPred_noCollab = log(1-(b + Xs_noCollab*beta'));

    log_likelihood = -(sum(linPred_collab) + sum(linPred_noCollab));

end
