clear all
close all
clc

% Load the datasets
gdpData = readtable('GDP.csv');
unemploymentData = readtable('Unemployment.csv');
inflationData = readtable('Inflation.csv');
resData = readtable('RES.csv');

varNames = string(table2cell(gdpData(1, :)));  % Extract the first row as a cell array and convert to string
varNames = cellstr(varNames);
gdpData.Properties.VariableNames = varNames;  % Set the variable names
unemploymentData.Properties.VariableNames = varNames;  % Set the variable names
inflationData.Properties.VariableNames = varNames;  % Set the variable names
resData.Properties.VariableNames = varNames;  % Set the variable names

% Define the desired country code
desiredCountryCode = 'POL'; % Replace with the desired country code

% Find the rows that correspond to the desired country in each dataset
countryRowsGDP = strcmp(gdpData.("Country Code"), desiredCountryCode);
rowIndexGDP = find(countryRowsGDP);

countryRowsUnemployment = strcmp(unemploymentData.("Country Code"), desiredCountryCode);
rowIndexUn = find(countryRowsUnemployment);

countryRowsInflation = strcmp(inflationData.("Country Code"), desiredCountryCode);
rowIndexIn = find(countryRowsInflation);

countryRowsRES = strcmp(resData.("Country Code"), desiredCountryCode);
rowIndexRES = find(countryRowsRES);

startColumnIndex = '1960';
endColumnIndex = '2022'; 

% Find the column indices based on the start and end column names
startColIdx = find(strcmp(gdpData.Properties.VariableNames, startColumnIndex));
endColIdx = find(strcmp(gdpData.Properties.VariableNames, endColumnIndex));

% Extract the rows for the desired country from each dataset
rowGDPData = gdpData(rowIndexGDP, startColIdx:endColIdx);
rowUnData = unemploymentData(rowIndexUn, startColIdx:endColIdx);
rowInData = inflationData(rowIndexIn, startColIdx:endColIdx);
rowRESData = resData(rowIndexRES, startColIdx:endColIdx);

% Convert the necessary variables to double
y = table2array(rowGDPData)';
x(:,1) = table2array(rowUnData)';
x(:,2) = table2array(rowInData)';
x(:,3) = table2array(rowRESData)';

% N = length(y);
% K = size(x,2);

% Remove rows with NaN values
nanRows = any(isnan(y), 2) | any(isnan(x), 2);
y(nanRows,:) = [];
x(nanRows,:) = [];
GDP0 = y;

N = size(y, 1);
K = size(x, 2);

%% Estimate the Smooth Transition Regression model with logistic smoothing function

s = x(:,1);

% grids
LAMBDA = logspace(-2,2,50); %0.01 = 10^(-2), 100 = 10^2, 50 values
C = linspace(quantile(s, 0.1), quantile(s, 0.9), 25); %25 values

% grid search
for i = 1:length(LAMBDA)
    for j = 1:length(C)
        lambda = LAMBDA(i);
        c = C(j);
        G = 1./(1+exp(-lambda*(s-c)));
        betas = OLS([x G.*x], G);
        params = [betas;lambda;c];
        LS(i,j) = loss2(y,x,s,params);
    end
end

[ii, jj] = find(LS==min(min(LS)));
lambda = LAMBDA(ii)
c = C(jj)
G = 1./(1+exp(-lambda*(s-c)));
betas = OLS([x G.*x], G);
beta1 = betas(1:K);
beta2 = betas(K+1:end);
params = [beta1;beta2;lambda;c];

% function handle
fun = @(param) loss2(y,x,s,param);
L = fun(params)

lb = [-inf(1,2*K), 0.01, quantile(s,0.1)];
ub = [inf(1,2*K), 100, quantile(s, 0.9)];
[param2, L1] = fmincon(fun, params, [], [], [], [], lb, ub);

% divide params
beta1 = param2(1:K);
beta2 = param2(K+1:2*K);
lam = param2(2*K+1);
c = param2(2*K+2);

GDP1 = x * beta1 + G .* x * beta2;

%% Estimate the Smooth Transition Regression model with exponential smoothing function

s = x(:,1);

% grids
LAMBDA = logspace(-2,2,50); %0.01 = 10^(-2), 100 = 10^2, 50 values
C = linspace(quantile(s, 0.1), quantile(s, 0.9), 25); %25 values

% grid search
for i = 1:length(LAMBDA)
    for j = 1:length(C)
        lambda = LAMBDA(i);
        c = C(j);
        G = 1./(1+exp(-lambda*(s-c)));
        betas = OLS([x G.*x], G);
        params = [betas;lambda;c];
        LS(i,j) = loss(x,y,s,params);
    end
end

[ii, jj] = find(LS==min(min(LS)));
lambda = LAMBDA(ii)
c = C(jj)
G2 = 1-exp(-lambda*(s-c).^2);
betas = OLS([x G2.*x], G2);
beta1 = betas(1:K);
beta2 = betas(K+1:end);
params = [beta1;beta2;lambda;c];

% function handle
fun = @(param) loss(x,y,s,param);
L = fun(params)

lb = [-inf(1,2*K), 0.01, quantile(s,0.1)];
ub = [inf(1,2*K), 100, quantile(s, 0.9)];
[param2, L2] = fmincon(fun, params, [], [], [], [], lb, ub);

% divide params
beta1 = param2(1:K);
beta2 = param2(K+1:2*K);
lam = param2(2*K+1);
c = param2(2*K+2);

GDP2 = x * beta1 + G2 .* x * beta2;

%% Linear Regression model

[betas3 L3] = OLS(x,G);
GDP3 = x * beta1;

%% Plotting

Gpl = table2array(rowGDPData)';
xpl(:,1) = table2array(rowUnData)';
xpl(:,2) = table2array(rowInData)';
xpl(:,3) = table2array(rowRESData)';

varNames = rowGDPData.Properties.VariableNames;
varNamesRow = string(varNames);

% Remove rows with NaN values
nanRows = any(isnan(Gpl), 2) | any(isnan(xpl), 2);
Gpl(nanRows,:) = [];
xpl(nanRows,:) = [];
varNamesRow(nanRows) = [];

years = str2double(varNamesRow');

% rowGDPData = rowGDPData(:, ~any(isnan(table2array(rowGDPData)), 1));
% varNames = rowGDPData.Properties.VariableNames;
% years = str2double(varNames);

figure;
hold on;
plot(years, GDP0, 'k', 'LineWidth', 2);
plot(years, GDP1, 'r--', 'LineWidth', 1.5);
plot(years, GDP2, 'g--', 'LineWidth', 1.5);
plot(years, GDP3, 'b--', 'LineWidth', 1.5);

% Customize the plot
legend('Real GDP', 'GDP (Logistic STR)', 'GDP (Exponential STR)', 'GDP (Linear Regression)', 'Location', 'best');
xlabel('Year');
ylabel('GDP');
xticks(years);
xticklabels(years);
xtickangle(45);

hold off;

%% BIC

N = size(G, 1);
K = size(x, 2);

BIC_logistic = bic2(L1, N, K);

BIC_exponential = bic2(L2, N, K);

BIC_linear = bic2(L3, N, K);

modelNames = {'Logistic STR', 'Exponential STR', 'Linear Regression'};
BIC_values = [BIC_logistic; BIC_exponential; BIC_linear];
BIC_table = table(modelNames', BIC_values, 'VariableNames', {'Model', 'BIC'});

disp(BIC_table)