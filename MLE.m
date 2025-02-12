

%% Read in the files  

filePath = '/Users/emma/Documents/GitHub/higher-order-collaborations/data/MCL_2015_matrix.csv';
T_MCL1 = readtable(filePath);

filePath = '/Users/emma/Documents/GitHub/higher-order-collaborations/data/TDA_2015_matrix.csv';
T_TDA1 = readtable(filePath);

filePath = '/Users/emma/Documents/GitHub/higher-order-collaborations/data/AES_2017_matrix.csv';
T_AES1 = readtable(filePath);

filePath = '/Users/emma/Documents/GitHub/higher-order-collaborations/data/CMC_2018_matrix.csv';
T_CMC1 = readtable(filePath);

allData = [T_MCL1; T_TDA1; T_AES1; T_CMC1];

%% Testing interaction for collaborators and non-collaborators at each conference

% Define conferences
conferences = {'MCL', 'TDA', 'AES', 'CMC'};

% Initialize table
T_summary = table();

% Loop over each conference
for i = 1:length(conferences)
    conf = conferences{i};
    T_data = eval(['T_', conf, '1']);  % Dynamically access table
    
    % Identify collaboration and non-collaboration rows
    collab_rows = T_data.collaborated == 1;
    no_collab_rows = T_data.collaborated == 0;

    num_collab = sum(collab_rows);
    num_no_collab = sum(no_collab_rows);

    % Compute mean values for sync1 to async3
    mean_sync1_collab = sum(T_data{collab_rows, "sync1"}) / num_collab;
    mean_sync1_no_collab = sum(T_data{no_collab_rows, "sync1"}) / num_no_collab;
    
    mean_async1_collab = sum(T_data{collab_rows, "async1"}) / num_collab;
    mean_async1_no_collab = sum(T_data{no_collab_rows, "async1"}) / num_no_collab;
    
    mean_async2_collab = sum(T_data{collab_rows, "async2"}) / num_collab;
    mean_async2_no_collab = sum(T_data{no_collab_rows, "async2"}) / num_no_collab;
    
    mean_async3_collab = sum(T_data{collab_rows, "async3"}) / num_collab;
    mean_async3_no_collab = sum(T_data{no_collab_rows, "async3"}) / num_no_collab;

    % Append results to table
    T_summary = [T_summary; 
        table({conf}, 1, mean_sync1_collab, mean_async1_collab, mean_async2_collab, mean_async3_collab, ...
              'VariableNames', {'Conference', 'Collab', 'sync1', 'async1', 'async2', 'async3'});
        table({conf}, 0, mean_sync1_no_collab, mean_async1_no_collab, mean_async2_no_collab, mean_async3_no_collab, ...
              'VariableNames', {'Conference', 'Collab', 'sync1', 'async1', 'async2', 'async3'})];
end

%% Aggregate across all tables 

% Initialize aggregated table
T_agg_summary = table();

% Identify collaboration and non-collaboration rows across all conferences
all_collab_rows = T_summary.Collab == 1;
all_no_collab_rows = T_summary.Collab == 0;

num_collab = sum(all_collab_rows);
num_no_collab = sum(all_no_collab_rows);

% Compute mean values for sync1 to async3 across all conferences
mean_sync1_collab = mean(T_summary{all_collab_rows, "sync1"});
mean_sync1_no_collab = mean(T_summary{all_no_collab_rows, "sync1"});

mean_async1_collab = mean(T_summary{all_collab_rows, "async1"});
mean_async1_no_collab = mean(T_summary{all_no_collab_rows, "async1"});

mean_async2_collab = mean(T_summary{all_collab_rows, "async2"});
mean_async2_no_collab = mean(T_summary{all_no_collab_rows, "async2"});

mean_async3_collab = mean(T_summary{all_collab_rows, "async3"});
mean_async3_no_collab = mean(T_summary{all_no_collab_rows, "async3"});

% Append results to aggregated table
T_agg_summary = [T_agg_summary; 
    table({"All Conferences"}, 1, mean_sync1_collab, mean_async1_collab, mean_async2_collab, mean_async3_collab, ...
          'VariableNames', {'Conference', 'Collab', 'sync1', 'async1', 'async2', 'async3'});
    table({"All Conferences"}, 0, mean_sync1_no_collab, mean_async1_no_collab, mean_async2_no_collab, mean_async3_no_collab, ...
          'VariableNames', {'Conference', 'Collab', 'sync1', 'async1', 'async2', 'async3'})];

% Display aggregated table
disp(T_agg_summary);


%%

% Extract unique conferences
confNames = unique(T_summary.Conference, 'stable');

% Extract values for plotting
VarNames = {'sync1', 'async1', 'async2', 'async3'};
numVars = length(VarNames);

% Reshape data for grouped bar chart
collab_vals = zeros(length(confNames), numVars);
no_collab_vals = zeros(length(confNames), numVars);

for i = 1:length(confNames)
    collab_vals(i, :) = T_summary{strcmp(T_summary.Conference, confNames{i}) & T_summary.Collab == 1, 3:6};
    no_collab_vals(i, :) = T_summary{strcmp(T_summary.Conference, confNames{i}) & T_summary.Collab == 0, 3:6};
end

% Create bar plots
figure;
for i = 1:numVars
    subplot(4,1,i);
    b = bar([collab_vals(:,i), no_collab_vals(:,i)], 'grouped');
    b(1).FaceColor = [0 0.4470 0.7410]; % MATLAB default blue (collab)
    b(2).FaceColor = [0.8500 0.3250 0.0980]; % MATLAB default red (no_collab)
    
    set(gca, 'XTickLabel', confNames);
    ylabel(VarNames{i});
    
    if i == 1
        legend({'Collaborators', 'Non-collaborators'}, 'Location', 'Best');
    end
end
xlabel('Conference');

%% 

% Extract values for plotting
VarNames = {'sync1', 'async1', 'async2', 'async3'};
numVars = length(VarNames);

collab_vals = 60*T_agg_summary{T_agg_summary.Collab == 1, 3:6}; % Mean values for collabs
no_collab_vals = 60*T_agg_summary{T_agg_summary.Collab == 0, 3:6}; % Mean values for non-collabs

% Create aggregated bar plot
figure;
barData = [collab_vals; no_collab_vals]';  % Transpose for proper bar grouping

b = bar(barData, 'grouped'); % Grouped bars (Collaborators vs Non-Collaborators)

% Fix bar colors
b(1).FaceColor = [0 0.4470 0.7410]; % MATLAB default blue (collab)
b(2).FaceColor = [0.8500 0.3250 0.0980]; % MATLAB default red (non-collab)

% FormattingInterIn
set(gca, 'FontSize', 18); % Set x-axis labels as variable names
ylabel('Interaction [min]','FontSize',18);
xticklabels({'Synchronous interaction','Asynchronous clique','Asynchronous wedge','Asynchronous link'});
legend({'Collaborators', 'Non-Collaborators'}, 'Location', 'Best');
grid off;


%%
rows1 = allData.collaborated == 1;
rows0 = allData.collaborated == 0;

sumcollaboratedeq1 = sum( sum( allData{rows1,1:4}, 2 ) );
sumcollaboratedeq0 = sum( sum( allData{rows0,1:4}, 2 ) );

count1 = sum(rows1);
count0 = sum(rows0);

avgcollaboratedeq1 = sumcollaboratedeq1 / count1;
avgcollaboratedeq0 = sumcollaboratedeq0 / count0;

%% 

rows1 = allData.collaborated == 1;
rows0 = allData.collaborated == 0;

data1 = sum(allData{rows1,1:4},2);
data0 = sum(allData{rows0,1:4},2);

statfun = @(x) mean(x);

nBoot = 1000;

[CI1, bootstats1] = bootci(nBoot,{statfun,data1},'type','bca');
[CI0, bootstats0] = bootci(nBoot,{statfun,data0},'type','bca');

avgcollaboratedeq1 = mean(data1);
avgcollaboratedeq0 = mean(data0);

%% 
avgVals = [avgcollaboratedeq1; avgcollaboratedeq0];
errLow = [avgcollaboratedeq1 - CI1(1); avgcollaboratedeq0 - CI0(1)];
errHigh = [CI1(2) - avgcollaboratedeq1; CI0(2) - avgcollaboratedeq0];

figure;
b = bar(1:2, avgVals, 'FaceColor','flat');
b.CData(1,:) = [0 0.4470 0.7410];   % Default MATLAB blue
b.CData(2,:) = [0.8500 0.3250 0.0980]; % Default MATLAB red
hold on;

errorbar(1:2, avgVals, errLow, errHigh, 'k', 'LineStyle', 'none');
set(gca, 'XTick', 1:2, 'XTickLabel', {'Collaborators','Non-collaborators'});
ylabel('Interaction [min]');

%% MLE with sync1

I_collab = allData.sync1(allData.collaborated == 1) ;
I_no_collab = allData.sync1(allData.collaborated == 0);

options = optimset('PlotFcns',@optimplotfval,'MaxFunEvals',1e4);
a = 0.0001;
b = 0.0001;
x0 = [a b];
LLwrapper = @(param) ( LL_2param(param(1),param(2),I_collab,I_no_collab)) ;
[x,fval] = fminsearch(LLwrapper,x0,options);
logL = fval;
AIC_sync = 2*2 +2*logL;

%% MLE with async1

collab_values = allData.async3(allData.collaborated == 1) ;
no_collab_values = allData.async3(allData.collaborated == 0);

options = optimset('PlotFcns',@optimplotfval,'MaxFunEvals',1e4);
a = 0.001;
b = 0.001;
x0 = [a b];
LLwrapper = @(param) ( LL_2param(param(1),param(2),collab_values,no_collab_values)) ;
[x,fval] = fminsearch(LLwrapper,x0,options);
logL = fval;
AIC = 2*2 + 2*logL

%% MLE for model where it is a constant probability between 0 and 1. 
options = optimset('PlotFcns',@optimplotfval,'MaxFunEvals',1e20);
p0 = 0.001;
LLwrapper = @(param) ( LL1(param(1),I_collab,I_no_collab)) ;
[x,fval] = fminsearch(LLwrapper,p0,options);

logL = fval;
AIC_const = 2*1 + 2*logL

%% Logistic regression:

T = allData;
X = T{:,1:4};  % sync1 to async3
y = T.collaborated;    % collaborated
mdlLogistic = fitglm(X,y,'Distribution','binomial'); 

%% LASSO logistic regression:
[B,FitInfo] = lassoglm(X,y,'binomial','CV',10); 
idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
bestLambda = FitInfo.Lambda(idxLambdaMinDeviance);
coefLasso = [FitInfo.Intercept(idxLambdaMinDeviance); B(:,idxLambdaMinDeviance)];

%% Testing the MLEs for all possible combinations 

paramNames = {'sync1','async1','async2','async3'};

% This creates all subsets (including the empty one).
allSubsets = {};
for k = 0:length(paramNames)
    % nchoosek(1:4, k) gives all ways to pick k elements out of 4
    combos = nchoosek(1:4, k);
    for row = 1:size(combos,1)
        allSubsets{end+1} = combos(row,:);
    end
end

% Suppose our entire dataset is in `allData`,
% with columns sync1, async1, async2, async3, collaborated
Xfull = [allData.sync1, allData.async1, allData.async2, allData.async3];
y = allData.collaborated;  % 0 or 1

results = table;  % to store results

modelIdx = 0;
for s = 1:length(allSubsets)
%for s = 1:1
    idxSubset = allSubsets{s};  % e.g. [1 3] means sync1 & async2, etc.

    % Number of parameters: 1 intercept + length(idxSubset) betas
    k = 1 + length(idxSubset);

    % Provide some initial guess for fminsearch
    % e.g., small random or zeros
    initialGuess = 0.0001*ones(1, k);

    lbnd = zeros(1,k); %lower bounds
    ubnd = 10*ones(1,k);

    % Create the wrapper function for fminsearch:
    LLwrapper = @(param) logisticNegLogLike(param, allData, idxSubset);

    % Fit by searching for param vector that minimizes negative log-likelihood
    [bestParam, fval] = fminsearchbnd(LLwrapper, initialGuess,lbnd,ubnd);

    % fval is the negative log-likelihood at the optimum
    negLogL = fval;

    % AIC = 2*k + 2*(-logL) 
    % Here negLogL = -logL, so AIC = 2*k + 2*negLogL
    AIC = 2*k + 2*negLogL;

    % Store results in a table or structure
    modelIdx = modelIdx + 1;
    results.modelIdx(modelIdx) = modelIdx;
    results.subset{modelIdx}   = idxSubset;      % which columns
    results.k(modelIdx)        = k;              % number of parameters
    results.negLogL(modelIdx)  = negLogL;
    results.AIC(modelIdx)      = AIC;
    results.params{modelIdx}   = bestParam;      % fitted parameter values
end

% Now you can examine `results` and find the best (lowest) AIC
[bestAIC, iBest] = min(results.AIC);
disp('Best AIC model:')
disp(results(iBest,:))

%% 
mdl = stepwiseglm( ...
    allData_new, ...
    'ResponseVar','collaborated', ...     % The dependent variable
    'PredictorVars',{'sync1','async1','async2','async3'}, ...
    'Distribution','binomial', ...        % Logistic regression
    'Link','logit', ...          
    'Upper','linear', ...                 % Full model: all main effects
    'Lower','constant', ...               % Minimal model: only intercept
    'Criterion','AIC', ...                % Use AIC for entering terms
    'Verbose',1);                         % Show step details

disp(mdl)              % Summaries of the final model
mdl.Coefficients       % Table of estimated coefficients
mdl.ModelCriterion     % AIC, BIC, etc.

%% 

%                 sync1     async1    async2    async3
oddsPerMinute = [1.3615, 1.0896, 1.2123, 1.3024];
labels        = {'sync1','async1','async2','async3'};

% Create a vector of times in minutes
t = 0:60;   % 61 points, from 0 to 60 inclusive

figure('Color','w');
hold on;

for i = 1:length(oddsPerMinute)
    % Compute the odds ratio at each minute: OR(t) = (OR_perMinute)^t
    or_t = oddsPerMinute(i).^t;
    
    % Plot a line for this interaction type
    plot(t, or_t, 'LineWidth',2, 'DisplayName',labels{i});
end

% Optional: if the curves become very large, consider using a log scale:
% set(gca,'YScale','log');

xlabel('Minutes');
ylabel('Odds Ratio');
title('Odds Ratio vs. Minutes');
legend('Location','Best');
set(gca,'FontSize',14);
grid on;
hold off;
