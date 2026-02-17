%% BLOCKS 
% signal used for comparison

%% Matrix Product case - many (M) number of simulations - using SNR
% Comparison of single base wavelet matrix with Matrix Products 
%% First: Two-Matrix Product W1W2 comparison to Single base
%% Second: Similarity Transform W1'W2W1 comparison to W1W2, Single base
%% Third: Kronecker Product Transform comparison to single base

% Three methods of Universal Thresholding, BAMS, and ABE are applied to the
% transformed signal  


close all
clear; clc;

% User Defined Parameters  
M = 200;  % number of simulations
SNR = 5;  % can change to 3, 5, or 7 


% Storing MSE values

% Kronecker Product for W0
mse_kron_vec = zeros(M, 1); 

% Matrix Products for W1'W2W1 
mse_W1W2W1_vec = zeros(M, 1); 

% Matrix Products for 2 filters 
mse_W1W2_vec = zeros(M, 1); 
 

% Single base
mse_W1_vec = zeros(M, 1);
mse_W2_vec = zeros(M, 1);


% BAMS vector initialization
mse_W1W2_BAMS_vec    = zeros(M,1);
mse_121_BAMS_vec     = zeros(M,1);
mse_W1_BAMS_vec      = zeros(M,1);
mse_W2_BAMS_vec      = zeros(M,1);
mse_kron_BAMS_vec    = zeros(M,1);


% ABE vector initialization
mse_W1W2_ABE_vec    = zeros(M,1);
mse_121_ABE_vec     = zeros(M,1);
mse_W1_ABE_vec      = zeros(M,1);
mse_W2_ABE_vec      = zeros(M,1);
mse_kron_ABE_vec    = zeros(M,1);

sigma2 = 1;   % because noise ~ N(0,1) and W is orthogonal

n=1024;
n1 = n;

% Blocks
m=n1;
t = (1:m) ./m;
pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
hgt = [4 (-5) 3 (-4) 5 (-4.2) 2.1 4.3  (-3.1) 2.1 (-4.2)];
sig = zeros(size(t));
for j=1:length(pos)
    sig = sig + (1 + sign(t-pos(j))).*(hgt(j)/2) ;
end

s_orig=sig;

figure(1)
plot(t, s_orig)

% Filter vectors 
hfilt_haar = [sqrt(2)/2, sqrt(2)/2]; % Haar filter


nl=3; % if using nl = 10 we run into problem with wavmat (W1 and W2 are not orthogonal!)       
W1= Wavmat(hfilt_haar, n, nl, 0); % Wavelet matrix 
W2= Wavmat(hfilt_haar, n, nl, 0); % Wavelet matrix 
 
% Add Gaussian noise with fixed SNR (user can set at beginning of file)
signal_var = var(s_orig);
signal_std = sqrt(signal_var);
s_orig= s_orig/signal_std * sqrt(SNR); % signal standardized to have unit variance, then scaled by sqrt(SNR)
     
for m = 1:M 
     s = s_orig + randn(size(s_orig)); % adds unit variance Gaussian noise
     
     %%  W1'W2W1 Product Wavelet Shrinkage 
     W121 = W1' * W2 * W1;
     wd121 = W121 * s'; % apply wavelet to noisy signal
     wd_121_raw = wd121;          % keep raw copy

     % Universal thresholding, noise ~ Normal(0,1)
     lambda = sqrt(2 * log(n)); 
    
     wd121(abs(wd121) < lambda) = 0;  % Hard thresholding
    
     sr121 = W121' * wd121; % reconstruct denoised signal

     % BAMS (bayesrule) shrinkage
     % --- uncomment and run this line for the levelwise BAMS---
     %wd_121_b = bamsShrinkCoeffs_Wavmat(wd_121_raw, n, nl);
     wd_121_b = bamsShrinkGlobal(wd_121_raw); % comment this global method line if running levelwise line above
     sr121_b = W121' * wd_121_b;

     % ABE shrinkage
     % sigma2 = 1 because noise ~ N(0,1) and W is orthogonal
     wd121_abe = abeShrink(wd_121_raw, sigma2);
     sr121_abe = W121' * wd121_abe;

     mse_W1W2W1_vec(m) = mean((s_orig - sr121').^2); 
     mse_121_BAMS_vec(m) = mean((s_orig - sr121_b').^2);
     mse_121_ABE_vec(m) = mean((s_orig - sr121_abe').^2);


     
     %%  W1W2 Product Wavelet Shrinkage 
     W12 = W1 * W2;
     wd12 = W12 * s'; % apply wavelet to noisy signal
     wd12_raw = wd12;          % keep raw copy


     % Universal thresholding, noise ~ Normal(0,1)
     lambda = sqrt(2 * log(n)); 
    
     wd12(abs(wd12) < lambda) = 0;  % Hard thresholding
    
     sr12 = W12' * wd12; % reconstruct denoised signal

     % BAMS (bayesrule) shrinkage
     % --- uncomment and run this line for the levelwise BAMS---
     %wd12_b = bamsShrinkCoeffs_Wavmat(wd12_raw, n, nl);
     wd12_b = bamsShrinkGlobal(wd12_raw); % comment this global method line if running levelwise line above
     sr12_b = W12' * wd12_b;


     % ABE shrinkage
     % sigma2 = 1 because noise ~ N(0,1) and W is orthogonal
     wd12_abe = abeShrink(wd12_raw, sigma2);
     sr12_abe = W12' * wd12_abe;



     mse_W1W2_vec(m) = mean((s_orig - sr12').^2); 
     mse_W1W2_BAMS_vec(m)  = mean((s_orig - sr12_b').^2);
     mse_W1W2_ABE_vec(m) = mean((s_orig - sr12_abe').^2);


     %%  Single Base (W1 - Haar) Wavelet Shrinkage  

     wd_W1 = W1 * s'; % Apply single wavelet matrix W1 to noisy signal
     wd_W1_raw = wd_W1;          % keep raw copy

     % Hard thresolding
     wd_thresh_W1 = wd_W1;
     wd_thresh_W1(abs(wd_W1) < lambda) = 0; 
    
     sr_W1 = W1' * wd_thresh_W1; % reconstruct denoised signal
    
     % BAMS (bayesrule) shrinkage
     % --- uncomment and run this line for the levelwise BAMS---
     %wd_W1_b = bamsShrinkCoeffs_Wavmat(wd_W1_raw, n, nl);
     wd_W1_b = bamsShrinkGlobal(wd_W1_raw); % comment this global method line if running levelwise line above
     srW1_b = W1' * wd_W1_b;

     % ABE shrinkage
     % sigma2 = 1 because noise ~ N(0,1) and W is orthogonal
     wd_W1_abe = abeShrink(wd_W1_raw, sigma2);
     sr_W1_abe = W1' * wd_W1_abe;


     mse_W1_vec(m) = mean((s_orig' - sr_W1).^2); % mse for single wavelet W1
     mse_W1_BAMS_vec(m) = mean((s_orig - srW1_b').^2);
     mse_W1_ABE_vec(m) = mean((s_orig - sr_W1_abe').^2);


     %%  Single Base (W2 - Haar) Wavelet Shrinkage  

     wd_W2 = W2 * s'; % Apply single wavelet matrix W2 to noisy signal
     wd_W2_raw = wd_W2;          % keep raw copy


     wd_thresh_W2 = wd_W2;
     wd_thresh_W2(abs(wd_W2) < lambda) = 0; % Hard thresolding
    
     sr_W2 = W2' * wd_thresh_W2; % reconstruct denoised signal

     % BAMS (bayesrule) shrinkage
     % --- uncomment and run this line for the levelwise BAMS---
     %wd_W2_b = bamsShrinkCoeffs_Wavmat(wd_W2_raw, n, nl);
     wd_W2_b = bamsShrinkGlobal(wd_W2_raw); % comment this global method line if running levelwise line above

     srW2_b = W2' * wd_W2_b;
    
     mse_W2_vec(m) = mean((s_orig' - sr_W2).^2); % mse for single wavelet W2
     mse_W2_BAMS_vec(m) = mean((s_orig - srW2_b').^2);
     
    %%  Kronecker Product 
    % Exploration of kron(W0,W0) 
     W0 = Wavmat(hfilt_haar, 2^7, 4, 0); % Wavelet matrix 
     %W0_2 = Wavmat(hfilt1, 2^3, 1, 0); % with level 1, variability is smaller
     W0_2 = Wavmat([sqrt(2)/2 sqrt(2)/2], 2^3, 1, 0); % with level 1, variability is smaller
     
     WW = kron(W0,W0_2);
     wd_kron = WW * s'; % apply wavelet to noisy signal
     wd_kron_raw = wd_kron;          % keep raw copy

     wd_kron(abs(wd_kron) < lambda) = 0;  % Hard thresholding
        
     sr_kron = WW' * wd_kron; % reconstruct denoised signal

     % BAMS (bayesrule) shrinkage
     % --- uncomment and run this line for the levelwise BAMS---
     %wd_kron_b = bamsShrinkCoeffs_Wavmat(wd_kron_raw, n, nl);
     wd_kron_b = bamsShrinkGlobal(wd_kron_raw); % comment this global method line if running levelwise line above
     srkron_b = WW' * wd_kron_b;

     % ABE shrinkage
     % sigma2 = 1 because noise ~ N(0,1) and W is orthogonal
     wd_kron_abe = abeShrink(wd_kron_raw, sigma2);
     srkron_abe = WW' * wd_kron_abe;
    
     mse_kron_vec(m) = mean((s_orig' - sr_kron).^2); % mse for kronecker prod
     mse_kron_BAMS_vec(m) = mean((s_orig - srkron_b').^2);
     mse_kron_ABE_vec(m) = mean((s_orig - srkron_abe').^2);


 end

 

 % Calculate Average MSE for W1'W2W1 across M simulations 
 avg_mse_W1W2W1 = mean(mse_W1W2W1_vec);

 % Calculate Average MSE for W1W2 across M simulations 
 avg_mse_W1W2 = mean(mse_W1W2_vec);

 % Calculate Average MSE for single bases across M simulations
 avg_mse_W1 = mean(mse_W1_vec);
 avg_mse_W2 = mean(mse_W2_vec);

 % Calculate Average MSE for kron(W0, W0_2) across M simulations 
 avg_mse_kron = mean(mse_kron_vec);

 % Calculate Average MSEs for all BAMS shrinkage-applied transforms
 avg_mse_W1_BAMS = mean(mse_W1_BAMS_vec);
 avg_mse_W2_BAMS = mean(mse_W2_BAMS_vec);
 avg_mse_W1W2_BAMS = mean(mse_W1W2_BAMS_vec);
 avg_mse_kron_BAMS = mean(mse_kron_BAMS_vec);
 avg_mse_121_BAMS = mean(mse_121_BAMS_vec);

 % Calculate Average MSEs for all ABE shrinkage-applied transforms
 avg_mse_W1_ABE = mean(mse_W1_ABE_vec);
 %avg_mse_W2_BAMS = mean(mse_W2_BAMS_vec); % W1 = W2 is Haar for blocks 
 avg_mse_W1W2_ABE = mean(mse_W1W2_ABE_vec);
 avg_mse_kron_ABE = mean(mse_kron_ABE_vec);
 avg_mse_121_ABE = mean(mse_121_ABE_vec);

 fprintf( 'ABE THRESHOLDING------------------------------------------------------\n'); 
 fprintf('| %-22s | %-11s | %-16s |\n', 'Method', 'Average MSE', 'Variance of MSE');
 fprintf('|------------------------|-------------|------------------|\n');

 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1 Haar (ABE)', avg_mse_W1_ABE, var(mse_W1_ABE_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2 (ABE)', avg_mse_W1W2_ABE, var(mse_W1W2_ABE_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'Kron (ABE)', avg_mse_kron_ABE, var(mse_kron_ABE_vec));

 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2W1 (ABE)', avg_mse_121_ABE, var(mse_121_ABE_vec));
 fprintf( '-----------------------------------------------------------------------\n'); 
 fprintf('\n')

 fprintf( 'BAMS THRESHOLDING-----------------------------------------------------\n'); 
 fprintf('| %-22s | %-11s | %-16s |\n', 'Method', 'Average MSE', 'Variance of MSE');
 fprintf('|------------------------|-------------|------------------|\n');

 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1 Haar (BAMS)', avg_mse_W1_BAMS, var(mse_W1_BAMS_vec));
 %fprintf('| %-22s | %.4f      | %.4f           |\n', 'W2 Haar (BAMS)', avg_mse_W2_BAMS, var(mse_W2_BAMS_vec));

 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2 (BAMS)', avg_mse_W1W2_BAMS, var(mse_W1W2_BAMS_vec));

 fprintf('| %-22s | %.4f      | %.4f           |\n', 'Kron (BAMS)', avg_mse_kron_BAMS, var(mse_kron_BAMS_vec));

 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2W1 (BAMS)', avg_mse_121_BAMS, var(mse_121_BAMS_vec));
 fprintf( '-----------------------------------------------------------------------\n'); 
 fprintf('\n')

 fprintf( 'UNIVERSAL THRESHOLDING----------------------------------------\n'); 
 fprintf('| %-22s | %-11s | %-16s |\n', 'Method', 'Average MSE', 'Variance of MSE');
 fprintf('|------------------------|-------------|------------------|\n');
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'Kronecker Product', avg_mse_kron, var(mse_kron_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2W1 Product', avg_mse_W1W2W1, var(mse_W1W2W1_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2 Product', avg_mse_W1W2, var(mse_W1W2_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1 (Haar)', avg_mse_W1, var(mse_W1_vec));
 %fprintf('| %-22s | %.4f      | %.4f           |\n', 'W2 (Haar)', avg_mse_W2, var(mse_W2_vec));


 %% Boxplots of Average MSEs

 % Product Transform vs. Single basis - Universal Thresholding
 figure;
 boxplot([mse_W1W2_vec, mse_W1_vec, mse_W2_vec], ...
     {'W1W2', 'Haar', 'Haar'});
 ylabel('MSE');
 title('MSE Distributions for M = 200 Simulations');

 % All Composite Transforms vs Single Bases - Universal Thresholding
 figure;
 boxplot([mse_kron_vec, mse_W1W2W1_vec, mse_W1W2_vec, mse_W1_vec, mse_W2_vec], ...
     {'kron', 'W1W2W1', 'W1W2', 'Haar', 'Haar'});
 ylabel('MSE');
 title('MSE Distributions including Kron and W1W2W1');
 
 
 % All Composite Transforms vs Single Bases - BAMS
 figure;
 boxplot([mse_kron_BAMS_vec, mse_121_BAMS_vec, mse_W1W2_BAMS_vec, mse_W1_BAMS_vec, mse_W2_BAMS_vec], ...
     {'kron', 'W1W2W1', 'W1W2', 'Haar', 'Haar'});
 ylabel('MSE');
 title('MSE Distributions including Kron and W1W2W1, BAMS thresholding');


 % All Composite Transforms vs Single Bases - ABE
 figure;
 boxplot([mse_kron_ABE_vec, mse_121_ABE_vec, mse_W1W2_ABE_vec, mse_W1_ABE_vec], ...
     {'kron', 'W1W2W1', 'W1W2', 'Haar'});
 ylabel('MSE');
 title('MSE Distributions including Kron and W1W2W1, ABE thresholding');
