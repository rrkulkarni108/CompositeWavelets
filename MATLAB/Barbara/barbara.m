
% Comparison of Product Matrices W1W2 to Single Base
% Usage of Barbara image for comparison
% Displays 4 images: 
% (1) Original Barbara 512 x 512 Image
% (2) Barbara image with noise added (sigma can be modified by choice of user)
% (3) Denoised Barbara image using single wavelet basis (bases can be chosen by user- we choose most symmetric wavelets)
% (4) Denoise Barbara image using product wavelet bases (bases chosen by user)

% Also outputs to console (1) MSEs of a single simulation of both single and product bases and 
%                         (2) average MSEs and var(MSE) across 200 simulations with customizable sigma noise value

% We perform hard thresholding with the universal threshold from VisuShrink (Donoho and Johnstone, 1992a)


close all; clear; clc;
colormap gray

%  User defined parameters  
M = 200;  % number of simulations
sigma = 20; % interesting to see is at sigma = 500                                


% Storing MSE values

% Matrix Products for 2 filters 
mse_W1W2_vec = zeros(M, 1); 
 
% Single base
mse_prod = zeros(M, 1);
mse_single = zeros(M, 1);

% Display original 512 x 512 image
Afull = double(imread('barbara_gray.bmp'));   % 512x512
A = Afull;  

[n, m] = size(A); 
N = n*n;

figure(1); imagesc(A); axis image off; title(sprintf('Barbara %dx%d',n,n));


Y = A + sigma*randn(size(A)); % signal + noise
figure(2); imagesc(Y); axis image off; colormap gray;
title('With noise added');


% Filter Library
hfilt_daub3 = [0.332670552950165   0.806891509311400   0.459877502118228  ...
     -0.135011020010067  -0.085441273882042  0.035226291882017]; % Daub3
hfilt_symm4 = [-0.075765714789273, ...
     -0.029635527645999, ...
      0.497618667632015, ...
      0.803738751805916, ...
      0.297857795605605, ...
     -0.099219543576935, ...
     -0.012603967262038, ...
      0.032223100604043];
hfilt_haar = [sqrt(2)/2, sqrt(2)/2]; % Haar filter
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

hfilt_coif3 = [	.038580777748	-.126969125396	-.077161555496	...
				.607491641386	.745687558934	.226584265197	]; % coiflet 6 tap

% Single and product bases
W1 = Wavmat(hfilt_symm4, n, 2, 2);
W2 = Wavmat(hfilt_daub3, n, 2, 2);
W1W2 = W1*W2;                               % product base

tau = sigma * sqrt(2*log(N)); % universal threshold

% Single base 
C_single = W1 * Y * W1';
C_single_th = C_single .* (abs(C_single) > tau);
Yhat_single = W1' * C_single_th * W1;
mse_single = mean((Yhat_single(:) - A(:)).^2);

% Product Matrix W1W2
C_prod = W1W2 * Y * W1W2';
C_prod_th = C_prod .* (abs(C_prod) > tau);
Yhat_prod = W1W2' * C_prod_th * W1W2;
mse_prod = mean((Yhat_prod(:) - A(:)).^2);


%plot
figure(3); imagesc(Yhat_single); axis image off; colormap gray; title('Denoised: Single base W');
figure(4); imagesc(Yhat_prod);   axis image off; colormap gray; title('Denoised: Product base W1W2');

fprintf('\nMSE of one simulation (single base W):      %.3f\n', mse_single);
fprintf('MSE of one simulation (product base W1W2):  %.3f\n', mse_prod);


% perform 200 simulations and calculate average mse
for m = 1:M 
    Y = A + sigma*randn(size(A)); % signal + noise
    tau = sigma * sqrt(2*log(N)); % Universal thresholding, noise ~ Normal(0,1)
    

    %%  W1W2 Product Wavelet Shrinkage  (W1- Symm 4, W2- DAUB 3)    
    C_prod = W1W2 * Y * W1W2'; % apply wavelet to noisy signal
    C_prod_th = C_prod .* (abs(C_prod) > tau); % Universal thresholding, noise ~ Normal(0,1)
    Yhat_prod = W1W2' * C_prod_th * W1W2; % reconstruct denoised signal
    mse_prod(m) = mean((Yhat_prod(:) - A(:)).^2);

    %%  Single Base (W1 - Symm4) Wavelet Shrinkage  
    C_single = W1 * Y * W1';
    C_single_th = C_single .* (abs(C_single) > tau);
    Yhat_single = W1' * C_single_th * W1;
    mse_single(m) = mean((Yhat_single(:) - A(:)).^2);


end

% Calculate Average MSE for W1W2 across M simulations 
 avg_mse_W1W2 = mean(mse_prod);

% and for single base W1
 avg_mse_W1 = mean(mse_single);

 % Print results

fprintf('| %-22s | %-11s | %-16s |\n', 'Method', 'Average MSE', 'Variance of MSE');
fprintf('|------------------------|-------------|------------------|\n');
fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2 Product', avg_mse_W1W2, var(mse_prod));
fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1 (Symm4)', avg_mse_W1, var(mse_single));




