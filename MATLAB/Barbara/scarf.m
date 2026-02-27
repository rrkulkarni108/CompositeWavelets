% Comparison of Product Matrices W1W2 to Single Base
% Usage of Barbara image for comparison
% Displays 4 images: 
% (1) Original Barbara 512 x 512 Image
% (2) Barbara image with noise added (sigma can be modified by choice of user)
% (3) Denoised Barbara image using single wavelet basis (bases can be chosen by user- we choose most symmetric wavelets)
% (4) Denoise Barbara image using product wavelet bases (bases chosen by user)
%
% Also outputs to console 
% (1) MSEs of a single simulation of both single and product bases 
% (2) average MSEs and var(MSE) across M simulations with customizable sigma
%
% We perform hard thresholding with the universal threshold from
% VisuShrink (Donoho and Johnstone, 1992a)

close all; clear; clc;
colormap gray

%  user defined parameters  
M     = 200;   % number of simulations
sigma = 20;    % noise std

% storing MSE values for the whole image
mse_prod_full   = zeros(M, 1);
mse_single_full = zeros(M, 1);
mse_single_full2 = zeros(M, 1);

% storing MSE values - scarf patch
mse_prod_scarf   = zeros(M, 1);
mse_single_scarf = zeros(M, 1);
mse_single_scarf2 = zeros(M, 1);


% Load image
Afull = double(imread('barbara_gray.bmp'));   % 512x512
A = Afull;  

[n, m] = size(A); 
N = n*n;

figure(1); imagesc(A); axis image off; title(sprintf('Barbara %dx%d',n,n));



% tablecloth indices 
%row_idx = 260:410;   
%col_idx =  60:200;  

%row_idx = 260:400;   
%col_idx =  5:100;   



% scarf indices-  scarf patch 
row_idx = 175:320;  
col_idx =  320:465;

A_scarf = A(row_idx, col_idx);



% One example noisy image
Y = A + sigma*randn(size(A)); % signal + noise
figure(2); imagesc(Y); axis image off; colormap gray;
title('With noise added (full image)');
figure(11); imagesc(A_scarf); axis image off; colormap gray;
title('Original scarf patch');
figure(12); imagesc(Y(row_idx, col_idx)); axis image off;colormap gray; 
title('Noisy scarf patch');

% Filter Library
hfilt_daub3 = [0.332670552950165   0.806891509311400   0.459877502118228  ...
               -0.135011020010067  -0.085441273882042   0.035226291882017]; % Daub3

hfilt_symm4 = [-0.075765714789273, ...
               -0.029635527645999, ...
                0.497618667632015, ...
                0.803738751805916, ...
                0.297857795605605, ...
               -0.099219543576935, ...
               -0.012603967262038, ...
                0.032223100604043];

hfilt_haar = [sqrt(2)/2, sqrt(2)/2]; % Haar filter

hfilt_daub5 = [ ...
    0.0033357252854737713;
   -0.012580751999015526;
   -0.0062414902130117052;
    0.077571493840045713;
   -0.032244869585029520;
   -0.24229488706619015;
    0.13842814590110342;
    0.72430852843857444;
    0.60382926979718952;
    0.16010239797419290];  % Daub 5

hfilt_coif1 = [ ...
     0.038580777748  -0.126969125396  -0.077161555496 ...
     0.607491641386   0.745687558934   0.226584265197]; % Coiflet 6 tap
hfilt_coif3  =  [	-.003793512864	.007782596426	.023452696142	...
				-.065771911281	-.061123390003	.405176902410	...
				.793777222626	.428483476378	-.071799821619	...
				-.082301927106	.034555027573	.015880544864	...
				-.009007976137	-.002574517688	.001117518771	...
				.000466216960	-.000070983303	-.000034599773	];
hfilt_coif4 = [	.000892313668	-.001629492013	-.007346166328	...
				.016068943964	.026682300156	-.081266699680	...
				-.056077313316	.415308407030	.782238930920	...
				.434386056491	-.066627474263	-.096220442034	...
				.039334427123	.025082261845	-.015211731527	...
				-.005658286686	.003751436157	.001266561929	...
				-.000589020757	-.000259974552	.000062339034	...
				.000031229876	-.000003259680	-.000001784985	]; % Coif 24 tap


% Single and product bases
W1   = Wavmat(hfilt_symm4, n, 2, 2);
W2   = Wavmat(hfilt_coif3, n, 2, 2);
W1W2 = W1*W2;   % product base

tau = sigma * sqrt(2*log(N)); % universal threshold

% One simulation, full image, scarf MSE

C_single    = W1 * Y * W1';
C_single_th = C_single .* (abs(C_single) > tau);
Yhat_single = W1' * C_single_th * W1;

C_prod    = W1W2 * Y * W1W2';
C_prod_th = C_prod .* (abs(C_prod) > tau);
Yhat_prod = W1W2' * C_prod_th * W1W2;

% Full-image MSE
mse_single_full_one = mean((Yhat_single(:) - A(:)).^2);
mse_prod_full_one   = mean((Yhat_prod(:)   - A(:)).^2);

% Scarf-patch MSE
Yhat_single_scarf = Yhat_single(row_idx, col_idx);
Yhat_prod_scarf   = Yhat_prod(row_idx, col_idx);

mse_single_scarf_one = mean((Yhat_single_scarf(:) - A_scarf(:)).^2);
mse_prod_scarf_one   = mean((Yhat_prod_scarf(:)   - A_scarf(:)).^2);

% Plot denoised full images
figure(3); imagesc(Yhat_single); axis image off; colormap gray; 
title('Denoised full image: Single base W');
figure(4); imagesc(Yhat_prod);   axis image off; colormap gray; 
title('Denoised full image: Product base W1W2');

% Plot denoised scarf patches
figure(13); imagesc(Yhat_single_scarf); axis image off;colormap gray; 
title('Denoised scarf patch: Single base W');
figure(14); imagesc(Yhat_prod_scarf); axis image off; colormap gray; 
title('Denoised scarf patch: Product base W1W2');

fprintf('\n--- One simulation (full image) ---\n');
fprintf('MSE full (single base W):      %.3f\n', mse_single_full_one);
fprintf('MSE full (product base W1W2):  %.3f\n', mse_prod_full_one);

fprintf('\n--- One simulation (scarf patch only) ---\n');
fprintf('MSE scarf (single base W):     %.3f\n', mse_single_scarf_one);
fprintf('MSE scarf (product base W1W2): %.3f\n', mse_prod_scarf_one);

% ---------------
% Perform M simulations and calculate average MSE for both the full image
% and scarf/tablecloth subset, depending on what is commented above

for i = 1:M 
    Y = A + sigma*randn(size(A));  % signal + noise
    tau = sigma * sqrt(2*log(N));  % universal threshold
    
    % Product base
    C_prod    = W1W2 * Y * W1W2';
    C_prod_th = C_prod .* (abs(C_prod) > tau);
    Yhat_prod = W1W2' * C_prod_th * W1W2;
    
    % Single base - Symm4
    C_single    = W1 * Y * W1';
    C_single_th = C_single .* (abs(C_single) > tau);
    Yhat_single = W1' * C_single_th * W1;
    
    % Single base - Coif3
    C_single2    = W2 * Y * W2';
    C_single_th2 = C_single2 .* (abs(C_single2) > tau);
    Yhat_single2 = W2' * C_single_th2 * W2;

    % Full-image MSE
    mse_prod_full(i)   = mean((Yhat_prod(:)   - A(:)).^2);
    mse_single_full(i) = mean((Yhat_single(:) - A(:)).^2);
    mse_single_full2(i) = mean((Yhat_single2(:) - A(:)).^2);
    
    % Scarf MSE
    A_scarf            = A(row_idx, col_idx); % (constant, but cheap)
    Yhat_prod_scarf    = Yhat_prod(row_idx, col_idx);
    Yhat_single_scarf  = Yhat_single(row_idx, col_idx);
    Yhat_single_scarf2  = Yhat_single2(row_idx, col_idx);

    
    mse_prod_scarf(i)   = mean((Yhat_prod_scarf(:)   - A_scarf(:)).^2);
    mse_single_scarf(i) = mean((Yhat_single_scarf(:) - A_scarf(:)).^2);
    mse_single_scarf2(i) = mean((Yhat_single_scarf2(:) - A_scarf(:)).^2);
end

% averages
avg_mse_prod_full   = mean(mse_prod_full);
avg_mse_single_full = mean(mse_single_full);
avg_mse_single_full2 = mean(mse_single_full2);


avg_mse_prod_scarf   = mean(mse_prod_scarf);
avg_mse_single_scarf = mean(mse_single_scarf);
avg_mse_single_scarf2 = mean(mse_single_scarf2);

fprintf('\n| %-22s | %-11s | %-16s |\n', 'Method (full)', 'Average MSE', 'Variance of MSE');
fprintf('|------------------------|-------------|------------------|\n');
fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2 Product', avg_mse_prod_full,   var(mse_prod_full));
fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1 (Symm4)',   avg_mse_single_full, var(mse_single_full));
fprintf('| %-22s | %.4f      | %.4f           |\n', 'W2 (Coif3)',   avg_mse_single_full2, var(mse_single_full2));


fprintf('\n| %-22s | %-11s | %-16s |\n', 'Method (scarf)', 'Average MSE', 'Variance of MSE');
fprintf('|------------------------|-------------|------------------|\n');
fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2 Product', avg_mse_prod_scarf,   var(mse_prod_scarf));
fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1 (Symm4)',   avg_mse_single_scarf, var(mse_single_scarf));
fprintf('| %-22s | %.4f      | %.4f           |\n', 'W2 (Coif3)',   avg_mse_single_scarf2, var(mse_single_scarf2));

% --------------------------------------------------------------------
% Plot MSE over iterations for scarf patch
figure;
plot(1:M, mse_single_scarf, '-o', 'DisplayName','Single base Symm4 (scarf)');
hold on;
plot(1:M, mse_single_scarf2, '-o', 'DisplayName','Single base Coif3 (scarf)');
hold on;
plot(1:M, mse_prod_scarf,   '-s', 'DisplayName','Product base (scarf)');
xlabel('Iteration'); ylabel('MSE (scarf patch)');
title('MSE across iterations on Barbara scarf patch');
legend('Location','best'); grid on;

% 4 Panel Figure: Original, Noisy, Scarf 
figure(100); clf;

subplot(2,2,1);
imagesc(A_scarf); axis image off; colormap gray;
title('Original Scarf Patch');

subplot(2,2,2);
imagesc(Y(row_idx, col_idx)); axis image off;
title(sprintf('Noisy Scarf Patch (\\sigma = %d)', sigma));

subplot(2,2,3);
imagesc(Yhat_single_scarf); axis image off;
title('Denoised (Single Base W1)');

subplot(2,2,4);
imagesc(Yhat_prod_scarf); axis image off;
title('Denoised (Product Base W1W2)');

%sgtitle('Scarf Patch Denoising Comparison');

% Normalize each patch to [0,1] for clean saving
orig_scarf   = mat2gray(A_scarf);
noisy_scarf  = mat2gray(Y(row_idx, col_idx));
single_scarf = mat2gray(Yhat_single_scarf);
prod_scarf   = mat2gray(Yhat_prod_scarf);

% Save as 300 dpi PNGs
%imwrite(orig_scarf,   'scarf_original.png');
%imwrite(noisy_scarf,  'scarf_noisy.png');
%imwrite(single_scarf, 'scarf_single.png');
%imwrite(prod_scarf,   'scarf_product.png');
