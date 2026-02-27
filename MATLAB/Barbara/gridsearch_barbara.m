%% Grid search of all filter pairs to determine which yields lowest MSE
%Comparison of Product Matrices W1W2 to Single Base
% Usage of Barbara image for comparison
% For practical purposes one may choose a shorter filter of the same type 
% since results can be very similar (such as Coif3 over Coif4)

close all; clear; clc;
colormap gray

%  User defined parameters  
M     = 200;   % number of simulations
sigma = 20;    % noise std

% Display original 512 x 512 image
Afull = double(imread('barbara_gray.bmp'));   % 512x512
A = Afull;  

[n, m] = size(A); 
N = n*n;

figure(1); imagesc(A); axis image off; title(sprintf('Barbara %dx%d',n,n));

% Example noisy image for visualization
Y = A + sigma*randn(size(A)); % signal + noise
figure(2); imagesc(Y); axis image off; colormap gray;
title(sprintf('With noise added (\\sigma = %.1f)', sigma));

% Filter Library

hfilt_daub3 = [0.332670552950165   0.806891509311400   0.459877502118228 ...
               -0.135011020010067  -0.085441273882042   0.035226291882017];
hfilt_symm4 = [-0.075765714789273, -0.029635527645999,  0.497618667632015, ...
                0.803738751805916,  0.297857795605605, -0.099219543576935, ...
               -0.012603967262038,  0.032223100604043];
                           % Symm4

hfilt_haar = [sqrt(2)/2, sqrt(2)/2];                            % Haar

hfilt_daub5 = [
    0.0033357252854737713;
   -0.012580751999015526;
   -0.0062414902130117052;
    0.077571493840045713;
   -0.032244869585029520;
   -0.24229488706619015;
    0.13842814590110342;
    0.72430852843857444;
    0.60382926979718952;
    0.16010239797419290
];  % Daub 5
hfilt_coif1 = [	.038580777748	-.126969125396	-.077161555496	...
				.607491641386	.745687558934	.226584265197	]; % coiflet 6 tap
hfilt_coif2 = [ 0.01638733646343463   -0.04146493678721403  -0.06737255472451776 ...
          0.3861100668241113     0.8127236354502723    0.4170051844212678  ...
         -0.07648859907841388   -0.0594344186451854    0.02368017194655153 ...
          0.005611434819108795  -0.001823208870779224 -0.0007205494455409002 ]; % coif 12 tap
hfilt_coif3  =  [	-.003793512864	.007782596426	.023452696142	...
				-.065771911281	-.061123390003	.405176902410	...
				.793777222626	.428483476378	-.071799821619	...
				-.082301927106	.034555027573	.015880544864	...
				-.009007976137	-.002574517688	.001117518771	...
				.000466216960	-.000070983303	-.000034599773	]; % coif 18 tap

hfilt_coif4 = [	.000892313668	-.001629492013	-.007346166328	...
				.016068943964	.026682300156	-.081266699680	...
				-.056077313316	.415308407030	.782238930920	...
				.434386056491	-.066627474263	-.096220442034	...
				.039334427123	.025082261845	-.015211731527	...
				-.005658286686	.003751436157	.001266561929	...
				-.000589020757	-.000259974552	.000062339034	...
				.000031229876	-.000003259680	-.000001784985	]; %coif 24 tap

% Put filters into a candidate list
cands = { ...
    struct('name','Haar',   'h', hfilt_haar), ...
    struct('name','Symm4',  'h', hfilt_symm4), ...
    struct('name','Db3',    'h', hfilt_daub3), ...
    struct('name','Db5',    'h', hfilt_daub5), ...
    struct('name','Coif1',  'h', hfilt_coif1) ...
    struct('name','Coif3',  'h', hfilt_coif3) ...
    struct('name','Coif4',  'h', hfilt_coif4) ...

};

K = numel(cands);

% Build wavelet matrices for all candidates (single bases)
Wcell = cell(K,1);
for k = 1:K
    Wcell{k} = Wavmat(cands{k}.h, n, 2, 2);
end

% Define methods: all single bases + all products W_i W_j

methods = struct('name',{},'T',{});

% Single bases
for i = 1:K
    methods(end+1).name = sprintf('%s (single)', cands{i}.name);
    methods(end).T      = Wcell{i};
end

% Product bases W_i W_j 
for i = 1:K
    for j = 1:K
        methods(end+1).name = sprintf('%s ∘ %s', cands{i}.name, cands{j}.name);
        methods(end).T      = Wcell{i} * Wcell{j};
    end
end

num_methods = numel(methods);


% Monte Carlo simulation for each method

tau = sigma * sqrt(2*log(N));      % universal hard threshold
mse_mat = zeros(num_methods, M);   

fprintf('Running grid search over %d methods, M=%d, sigma=%.1f...\n', ...
        num_methods, M, sigma);

for mrep = 1:M
    Y = A + sigma*randn(size(A));  % noisy image (new noise each time)
    
    for p = 1:num_methods
        T = methods(p).T;
        
        % Wavelet transform
        C = T * Y * T';
        % Universal hard threshold
        C_th = C .* (abs(C) > tau);
        % Reconstruction
        Yhat = T' * C_th * T;
        
        % MSE
        mse_mat(p, mrep) = mean((Yhat(:) - A(:)).^2);
    end
end


% Summarize results

avg_mse = mean(mse_mat, 2);
var_mse = var(mse_mat, 0, 2);

% Build table and sort by Average MSE in ascending order
Method   = string({methods.name})';
AverageMSE = avg_mse;
VarMSE     = var_mse;

results_tbl = table(Method, AverageMSE, VarMSE);
results_tbl = sortrows(results_tbl, 'AverageMSE', 'ascend');

fprintf('\n| %-25s | %-12s | %-12s |\n', 'Method', 'Average MSE', 'Var(MSE)');
fprintf('|%s|\n', repmat('-',1,25+3+12+3+12+2));
for i = 1:height(results_tbl)
    fprintf('| %-25s | %10.4f | %10.4f |\n', ...
        results_tbl.Method(i), results_tbl.AverageMSE(i), results_tbl.VarMSE(i));
end

disp(' ');
disp('Ranked methods (best to worst by Average MSE):');
disp(results_tbl);


