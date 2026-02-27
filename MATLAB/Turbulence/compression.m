%% Turbulence compression: W1 (Symm4), W2 (Coif3), W1W2
% We make a comparison of energy concentration across transforms using top-1% and top-5% coefficient energy fractions.
close all; clear; clc;

% Load turbulence data (vertical velocity time series)
data = load('turbulence.txt');
x = data(:, 3);      % colulmn 3 has the vertical velocity (W)
x = x(:);
N = 2048; % only use part of the data since if we use the whole set then our computations would blow up
if length(x) < N
    error('Signal length (%d) is shorter than N = %d.', length(x), N);
end
x = x(1:N);

fprintf('Using N = %d samples from turbulence data (column 3).\n', N);


% Filters we will use
hfilt_symm4 = [-0.075765714789273, ...
               -0.029635527645999, ...
                0.497618667632015, ...
                0.803738751805916, ...
                0.297857795605605, ...
               -0.099219543576935, ...
               -0.012603967262038, ...
                0.032223100604043];


hfilt_coif3  =  [	-.003793512864	.007782596426	.023452696142	...
				-.065771911281	-.061123390003	.405176902410	...
				.793777222626	.428483476378	-.071799821619	...
				-.082301927106	.034555027573	.015880544864	...
				-.009007976137	-.002574517688	.001117518771	...
				.000466216960	-.000070983303	-.000034599773	];


% Single and product bases
% Using the same last two parameters as in Barbara - 2,2
W1 = Wavmat(hfilt_symm4, N, 2, 2);   % Symm4
W2 = Wavmat(hfilt_coif3,  N, 2, 2);   % Coif3

% Composite: W1W2 = W1 * W2  
W1W2 = W1 * W2;

c1  = W1   * x;   % Symm4, single base coefs
c2  = W2   * x;   % Coif3, single base coefs
c12 = W1W2 * x;   % Product W1W2 coefs



% compute epsilons

% W1 (Symm4)
eps1_1 = energyFractions(c1, 0.01);
eps5_1 = energyFractions(c1, 0.05);
epsW1_top1 = eps1_1.E_top / eps1_1.E_total;
epsW1_top5 = eps5_1.E_top / eps5_1.E_total;

% W2 (Coif3)
eps1_2 = energyFractions(c2, 0.01);
eps5_2 = energyFractions(c2, 0.05);
epsW2_top1 = eps1_2.E_top / eps1_2.E_total;
epsW2_top5 = eps5_2.E_top / eps5_2.E_total;

% W1W2 product
eps1_12 = energyFractions(c12, 0.01);
eps5_12 = energyFractions(c12, 0.05);
epsW1W2_top1 = eps1_12.E_top / eps1_12.E_total;
epsW1W2_top5 = eps5_12.E_top / eps5_12.E_total;


%% Display table 
fprintf('\nCompression of Turbulent Time Series (N = %d)\n', N);
fprintf('---------------------------------------------------------------\n');
fprintf('Transform     TotalEnergy       Top1%%%%Energy    Top5%%%%Energy\n');
fprintf('---------------------------------------------------------------\n');
fprintf('W1 (Symm4)   %12.4f       %10.4f        %10.4f\n', ...
        eps1_1.E_total,  epsW1_top1,   epsW1_top5);
fprintf('W2 (Coif3)   %12.4f       %10.4f        %10.4f\n', ...
        eps1_2.E_total,  epsW2_top1,   epsW2_top5);
fprintf('W1W2         %12.4f       %10.4f        %10.4f\n', ...
        eps1_12.E_total, epsW1W2_top1, epsW1W2_top5);

% W2W1 instead of W1W2:
W2W1 = W2 * W1;
c21 = W2W1 * x;
eps1_21 = energyFractions(c21, 0.01);
eps5_21 = energyFractions(c21, 0.05);
epsW2W1_top1 = eps1_21.E_top / eps1_21.E_total;
epsW2W1_top5 = eps5_21.E_top / eps5_21.E_total;
fprintf('W2W1         %12.4f       %10.4f        %10.4f\n', ...
         eps1_21.E_total, epsW2W1_top1, epsW2W1_top5);


% function that calculates the energy fractions
function stats = energyFractions(c, p)
    % c : coefficient vector
    % p : proportion 
    E_total = sum(c.^2);
    k = max(1, round(p * numel(c)));   % at least 1 coefficient
    s = sort(c.^2, 'descend');
    E_top = sum(s(1:k));

    stats.E_total = E_total;
    stats.E_top   = E_top;
end
