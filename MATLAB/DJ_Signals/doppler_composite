%% DOPPLER
% signal used for comparison

%% Matrix Product case - many (M) number of simulations - using SNR
% Comparison of single base wavelet matrix with Matrix Products :
%% First: Two-Matrix Product W1W2 comparison to Single base
%% Second: Similarity Transform W1'W2W1 comparison to W1W2, Single base
%% Third: Kronecker Product Transform comparison to single base

% Three methods of Universal Thresholding, BAMS, and ABE are applied to the
% transformed signal  

close all
clear; clc;

%  User Defined Parameters  
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

% sigma2 for ABE 
sigma2 = 1;   % because noise ~ N(0,1) and W is orthogonal

n1=512;

% Signal: Doppler
n = 1024;
t = linspace(0, 1, n);
doppler = sqrt(t .* (1 - t)) .* sin((2 * pi * 1.05) ./ (t + 0.05));
s_orig = doppler;  

% Filter vectors 
hfilt1 = [
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

 hfilt2=[-0.075765714789273, ...
     -0.029635527645999, ...
      0.497618667632015, ...
      0.803738751805916, ...
      0.297857795605605, ...
     -0.099219543576935, ...
     -0.012603967262038, ...
      0.032223100604043]; % Symmlet 4 for Doppler

 hfilt3 = [	 -0.003382415951003908, ...
        -0.000542132331797018, ...
        0.03169508781151886, ...
        0.00760748732494897,...
       -0.1432942383512576,...
       -0.06127335906765891,...
        0.4813596512592537,...
        0.7771857516996492, ...
        0.3644418948360139, ...
        -0.05194583810802026, ...
        -0.02721902991713553,...
        0.04913717967372511, ...
        0.003808752013880463, ...
        -0.01495225833706814, ...
        -0.0003029205147226741, ...
        0.001889950332768561 	]; % Symmlet 8



 nl=3; % if using nl = 10 we run into problem with wavmat (W1 and W2 are not orthogonal!)       
 W1= Wavmat(hfilt1, n, nl, 0); % Wavelet matrix 
 W2= Wavmat(hfilt2, n, nl, 0); % Wavelet matrix 
 
% Add Gaussian noise with fixed SNR (set at beginning of file)
signal_var = var(s_orig);
signal_std = sqrt(signal_var);
s_orig= s_orig/signal_std * sqrt(SNR); % signal standardized to have unit variance, then scaled by sqrt(SNR)
     
for m = 1:M 
     s = s_orig + randn(size(s_orig)); % adds unit variance Gaussian noise

     
     %%  W1'W2W1 Product Wavelet Shrinkage 
     W121 = W1' * W2 * W1;
     wd121 = W121 * s'; % apply wavelet to noisy signal

     wd_121_raw = wd121;

     lambda = sqrt(2 * log(n)); % Universal thresholding, noise ~ Normal(0,1)
    
     wd121(abs(wd121) < lambda) = 0;  % Hard thresholding
    
     sr121 = W121' * wd121; % reconstruct denoised signal

     % BAMS (global)
     wd_121_b = bamsShrinkGlobal(wd_121_raw);
     sr121_b  = W121' * wd_121_b;

     %  ABE 
     wd_121_abe = abeShrink(wd_121_raw, sigma2);
     sr121_abe  = W121' * wd_121_abe;

     %(s_orig-sr121') * (s_orig-sr121')'
     mse_W1W2W1_vec(m) = mean((s_orig - sr121').^2); 
     % store BAMS/ABE MSEs
     mse_121_BAMS_vec(m) = mean((s_orig - sr121_b').^2);
     mse_121_ABE_vec(m)  = mean((s_orig - sr121_abe').^2);

     
     %%  W1W2 Product Wavelet Shrinkage 
     W12 = W1 * W2;
     wd12 = W12 * s'; % apply wavelet to noisy signal
     % keep raw copy for BAMS/ABE 
     wd12_raw = wd12;

     lambda = sqrt(2 * log(n)); % Universal thresholding, noise ~ Normal(0,1)
    
     wd12(abs(wd12) < lambda) = 0;  % Hard thresholding
    
     sr12 = W12' * wd12; % reconstruct denoised signal

     %  BAMS (global) shrinkage
     wd12_b = bamsShrinkGlobal(wd12_raw);
     sr12_b = W12' * wd12_b;

     %  ABE shrinkage
     wd12_abe = abeShrink(wd12_raw, sigma2);
     sr12_abe = W12' * wd12_abe;

     %(s_orig-sr12') * (s_orig-sr12')'
     mse_W1W2_vec(m) = mean((s_orig - sr12').^2); 
     % store BAMS/ABE MSEs 
     mse_W1W2_BAMS_vec(m) = mean((s_orig - sr12_b').^2);
     mse_W1W2_ABE_vec(m)  = mean((s_orig - sr12_abe').^2);


     %%  Single Base (W1 - Daub 5) Wavelet Shrinkage  

     wd_W1 = W1 * s'; % Apply single wavelet matrix W1 to noisy signal
     %  keep raw copy for BAMS/ABE 
     wd_W1_raw = wd_W1;

     wd_thresh_W1 = wd_W1;
     wd_thresh_W1(abs(wd_W1) < lambda) = 0; % Hard thresolding
    
     sr_W1 = W1' * wd_thresh_W1; % reconstruct denoised signal
     
     %   BAMS (global) 
     wd_W1_b = bamsShrinkGlobal(wd_W1_raw);
     srW1_b  = W1' * wd_W1_b;

     %   ABE 
     wd_W1_abe = abeShrink(wd_W1_raw, sigma2);
     srW1_abe  = W1' * wd_W1_abe;

     mse_W1_vec(m) = mean((s_orig' - sr_W1).^2); % mse for single wavelet W1
     %   store BAMS/ABE MSEs 
     mse_W1_BAMS_vec(m) = mean((s_orig - srW1_b').^2);
     mse_W1_ABE_vec(m)  = mean((s_orig - srW1_abe').^2);
    

     %%  Single Base (W2 - Symm4) Wavelet Shrinkage  

     wd_W2 = W2 * s'; % Apply single wavelet matrix W2 to noisy signal
     %   keep raw copy for BAMS/ABE 
     wd_W2_raw = wd_W2;

     wd_thresh_W2 = wd_W2;
     wd_thresh_W2(abs(wd_W2) < lambda) = 0; % Hard thresolding
    
     sr_W2 = W2' * wd_thresh_W2; % reconstruct denoised signal
     
     %   BAMS (global) 
     wd_W2_b = bamsShrinkGlobal(wd_W2_raw);
     srW2_b  = W2' * wd_W2_b;

     %   ABE 
     wd_W2_abe = abeShrink(wd_W2_raw, sigma2);
     srW2_abe  = W2' * wd_W2_abe;

     mse_W2_vec(m) = mean((s_orig' - sr_W2).^2); % mse for single wavelet W2
     
     %  store BAMS/ABE MSEs 
     mse_W2_BAMS_vec(m) = mean((s_orig - srW2_b').^2);
     mse_W2_ABE_vec(m)  = mean((s_orig - srW2_abe').^2);
     
    %%  Kronecker Product 
    % Exploration of kron(W0,W0) 
     W0 = Wavmat(hfilt1, 2^7, 4, 0); % Wavelet matrix 
     %W0_2 = Wavmat(hfilt1, 2^3, 1, 0); % with level 1, variability is smaller
     W0_2 = Wavmat([sqrt(2)/2 sqrt(2)/2], 2^3, 1, 0); % with level 1, variability is smaller

     WW = kron(W0,W0_2);
     wd_kron = WW * s'; % apply wavelet to noisy signal
     wd_kron_raw = wd_kron; % keep raw copy for BAMS/ABE 

     wd_kron(abs(wd_kron) < lambda) = 0;  % Hard thresholding
        
     sr_kron = WW' * wd_kron; % reconstruct denoised signal
     
     %   BAMS (global) 
     wd_kron_b = bamsShrinkGlobal(wd_kron_raw);
     srkron_b  = WW' * wd_kron_b;

     %   ABE 
     wd_kron_abe = abeShrink(wd_kron_raw, sigma2);
     srkron_abe  = WW' * wd_kron_abe;

     mse_kron_vec(m) = mean((s_orig' - sr_kron).^2); % mse for kronecker prod
     %   store BAMS/ABE MSEs 
     mse_kron_BAMS_vec(m) = mean((s_orig - srkron_b').^2);
     mse_kron_ABE_vec(m)  = mean((s_orig - srkron_abe').^2);

 end

 % Calculate Average MSE for kron(W0, W0_2) across M simulations 
 avg_mse_kron = mean(mse_kron_vec);

 % Calculate Average MSE for W1'W2W1 across M simulations 
 avg_mse_W1W2W1 = mean(mse_W1W2W1_vec);

 % Calculate Average MSE for W1W2 across M simulations 
 avg_mse_W1W2 = mean(mse_W1W2_vec);

 % Calculate Average MSE for all single bases across M simulations
 avg_mse_W1 = mean(mse_W1_vec);
 avg_mse_W2 = mean(mse_W2_vec);

 %  Summary stats for BAMS 
 avg_mse_kron_BAMS = mean(mse_kron_BAMS_vec);
 avg_mse_121_BAMS  = mean(mse_121_BAMS_vec);
 avg_mse_W1W2_BAMS = mean(mse_W1W2_BAMS_vec);
 avg_mse_W1_BAMS   = mean(mse_W1_BAMS_vec);
 avg_mse_W2_BAMS   = mean(mse_W2_BAMS_vec);

 %  Summary stats for ABE 
 avg_mse_kron_ABE = mean(mse_kron_ABE_vec);
 avg_mse_121_ABE  = mean(mse_121_ABE_vec);
 avg_mse_W1W2_ABE = mean(mse_W1W2_ABE_vec);
 avg_mse_W1_ABE   = mean(mse_W1_ABE_vec);
 avg_mse_W2_ABE   = mean(mse_W2_ABE_vec);


 % Print results 

 fprintf('\nABE THRESHOLDING------------------------------------------------------\n');
 fprintf('| %-22s | %-11s | %-16s |\n', 'Method', 'Average MSE', 'Variance of MSE');
 fprintf('|------------------------|-------------|------------------|\n');
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'Kronecker Product (ABE)', avg_mse_kron_ABE, var(mse_kron_ABE_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2W1 Product (ABE)',    avg_mse_121_ABE,  var(mse_121_ABE_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2 Product (ABE)',      avg_mse_W1W2_ABE, var(mse_W1W2_ABE_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1 Daub 5 (ABE)',         avg_mse_W1_ABE,   var(mse_W1_ABE_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W2 Symm4 (ABE)',          avg_mse_W2_ABE,   var(mse_W2_ABE_vec));
 fprintf('-----------------------------------------------------------------------\n');

 fprintf('\nBAMS THRESHOLDING-----------------------------------------------------\n');
 fprintf('| %-22s | %-11s | %-16s |\n', 'Method', 'Average MSE', 'Variance of MSE');
 fprintf('|------------------------|-------------|------------------|\n');
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'Kronecker Product (BAMS)', avg_mse_kron_BAMS, var(mse_kron_BAMS_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2W1 Product (BAMS)',    avg_mse_121_BAMS,  var(mse_121_BAMS_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2 Product (BAMS)',      avg_mse_W1W2_BAMS, var(mse_W1W2_BAMS_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1 Daub 5 (BAMS)',         avg_mse_W1_BAMS,   var(mse_W1_BAMS_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W2 Symm4 (BAMS)',          avg_mse_W2_BAMS,   var(mse_W2_BAMS_vec));
 fprintf('-----------------------------------------------------------------------\n');

 
 fprintf('\nUNIVERSAL THRESHOLDING------------------------------------------------------\n');
 fprintf('| %-22s | %-11s | %-16s |\n', 'Method', 'Average MSE', 'Variance of MSE');
 fprintf('|------------------------|-------------|------------------|\n');
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'Kronecker Product', avg_mse_kron, var(mse_kron_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2W1 Product', avg_mse_W1W2W1, var(mse_W1W2W1_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2 Product', avg_mse_W1W2, var(mse_W1W2_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1 (Daubechies 5)', avg_mse_W1, var(mse_W1_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W2 (Symmlet 4)', avg_mse_W2, var(mse_W2_vec));

 
 %% Boxplots of Average MSEs

 % Product Transform vs. Single basis - Universal Thresholding
 figure;
 boxplot([mse_W1W2_vec, mse_W1_vec, mse_W2_vec], ...
     {'W1W2', 'Daubechies 5', 'Symmlet 4'});
 ylabel('MSE');
 title('MSE Distributions for M = 200 Simulations');


 % All Composite Transforms vs Single Bases - Universal Thresholding
 figure;
 boxplot([mse_kron_vec, mse_W1W2W1_vec, mse_W1W2_vec, mse_W1_vec, mse_W2_vec], ...
     {'kron', 'W1W2W1', 'W1W2', 'Daubechies 5', 'Symmlet 4'});
 ylabel('MSE');
 title('MSE Distributions including Kron and W1W2W1');

 % All Composite Transforms vs Single Bases - BAMS
 figure;
 boxplot([mse_kron_BAMS_vec, mse_121_BAMS_vec, mse_W1W2_BAMS_vec, mse_W1_BAMS_vec, mse_W2_BAMS_vec], ...
     {'kron', 'W1W2W1', 'W1W2', 'Daub5', 'Symm4'});
 ylabel('MSE');
 title('MSE Distributions including Kron and W1W2W1, BAMS thresholding');

 % All Composite Transforms vs Single Bases - ABE
 figure;
 boxplot([mse_kron_ABE_vec, mse_121_ABE_vec, mse_W1W2_ABE_vec, mse_W1_ABE_vec, mse_W2_ABE_vec], ...
     {'kron', 'W1W2W1', 'W1W2', 'Daub5', 'Symm4'});
 ylabel('MSE');
 title('MSE Distributions including Kron and W1W2W1, ABE thresholding');
