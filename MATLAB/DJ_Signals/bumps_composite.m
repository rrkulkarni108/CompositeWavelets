%% BUMPS 

%% Matrix Product case - many (M) number of simulations - using SNR
% Comparison of single base wavelet matrix with Matrix Products 
%% First: Two-Matrix Product W1W2 comparison to Single base
%% Second: Similarity Transform W1'W2W1 comparison to W1W2, Single base
%% Third: Kronecker Product 


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


n=1024;
n1 = n;

% Bumps
m=n1;
t = (1:m) ./m;
pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
hgt = [ 4  5   3   4  5  4.2 2.1 4.3  3.1 5.1 4.2];
wth = [.005 .005 .006 .01 .01 .03 .01 .01  .005 .008 .005];
sig = zeros(size(t));
for j =1:length(pos)
	sig = sig + hgt(j)./( 1 + abs((t - pos(j))./wth(j))).^4;
end 
s_orig=sig;

figure(1)
plot(t, s_orig)
 

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

hfilt_daub3 = [0.332670552950165   0.806891509311400   0.459877502118228  ...
     -0.135011020010067  -0.085441273882042  0.035226291882017]; % Daub3

hfilt_haar = [sqrt(2)/2, sqrt(2)/2]; % Haar filter


nl=2; % if using nl = 10 we run into problem with wavmat (W1 and W2 are not orthogonal!)       
W1= Wavmat(hfilt_haar, n, nl, 0); % Wavelet matrix 
W2= Wavmat(hfilt_daub3, n, nl, 0); % Wavelet matrix 
 
% Add Gaussian noise with fixed SNR (user can set at beginning of file)
signal_var = var(s_orig);
signal_std = sqrt(signal_var);
s_orig= s_orig/signal_std * sqrt(SNR); % signal standardized to have unit variance, then scaled by sqrt(SNR)
     

for m = 1:M 
     s = s_orig + randn(size(s_orig)); % adds unit variance Gaussian noise

     
     %%  W1'W2W1 Product Wavelet Shrinkage 
     W121 = W1 * W2 * W1;
     wd121 = W121 * s'; % apply wavelet to noisy signal

     lambda = sqrt(2 * log(n)); % Universal thresholding, noise ~ Normal(0,1)
    
     wd121(abs(wd121) < lambda) = 0;  % Hard thresholding
    
     sr121 = W121' * wd121; % reconstruct denoised signal

     mse_W1W2W1_vec(m) = mean((s_orig - sr121').^2); 

     
     %%  W1W2 Product Wavelet Shrinkage 
     W12 = W1 * W2;
     wd12 = W12 * s'; % apply wavelet to noisy signal

     lambda = sqrt(2 * log(n)); % Universal thresholding, noise ~ Normal(0,1)
    
     wd12(abs(wd12) < lambda) = 0;  % Hard thresholding
    
     sr12 = W12' * wd12; % reconstruct denoised signal

     mse_W1W2_vec(m) = mean((s_orig - sr12').^2); 


     %%  Single Base (W1 - Daub 3) Wavelet Shrinkage  

     wd_W1 = W1 * s'; % Apply single wavelet matrix W1 to noisy signal

     wd_thresh_W1 = wd_W1;
     wd_thresh_W1(abs(wd_W1) < lambda) = 0; % Hard thresolding
    
     sr_W1 = W1' * wd_thresh_W1; % reconstruct denoised signal
    
     mse_W1_vec(m) = mean((s_orig' - sr_W1).^2); % mse for single wavelet W1
    

     %%  Single Base (W2 - Daub 3) Wavelet Shrinkage  

     wd_W2 = W2 * s'; % Apply single wavelet matrix W2 to noisy signal

     wd_thresh_W2 = wd_W2;
     wd_thresh_W2(abs(wd_W2) < lambda) = 0; % Hard thresolding
    
     sr_W2 = W2' * wd_thresh_W2; % reconstruct denoised signal
    
     mse_W2_vec(m) = mean((s_orig' - sr_W2).^2); % mse for single wavelet W2
     

    %%  Kronecker Product 

    % Exploration of kron(W0,W0) 

     W0 = Wavmat(hfilt_daub3, 2^7, 4, 0); % Wavelet matrix 
     %W0_2 = Wavmat(hfilt1, 2^3, 1, 0); % with level 1, variability is smaller
     W0_2 = Wavmat([sqrt(2)/2 sqrt(2)/2], 2^3, 1, 0); % with level 1, variability is smaller

     WW = kron(W0,W0_2);
     wd_kron = WW * s'; % apply wavelet to noisy signal
     wd_kron(abs(wd_kron) < lambda) = 0;  % Hard thresholding
        
     sr_kron = WW' * wd_kron; % reconstruct denoised signal
    
     %(s_orig-sr_kron') * (s_orig-sr_kron')'
     mse_kron_vec(m) = mean((s_orig' - sr_kron).^2); % mse for kronecker prod


 end

 % Calculate Average MSE for kron(W0, W0_2) across M simulations 
 avg_mse_kron = mean(mse_kron_vec);

 % Calculate Average MSE for W1'W2W1 across M simulations 
 avg_mse_W1W2W1 = mean(mse_W1W2W1_vec);

 % Calculate Average MSE for W1W2 across M simulations 
 avg_mse_W1W2 = mean(mse_W1W2_vec);


 % and all single bases 
 avg_mse_W1 = mean(mse_W1_vec);
 avg_mse_W2 = mean(mse_W2_vec);

 % Print results

 fprintf('| %-22s | %-11s | %-16s |\n', 'Method', 'Average MSE', 'Variance of MSE');
 fprintf('|------------------------|-------------|------------------|\n');
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'Kronecker Product', avg_mse_kron, var(mse_kron_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2W1 Product', avg_mse_W1W2W1, var(mse_W1W2W1_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2 Product', avg_mse_W1W2, var(mse_W1W2_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1 (Daubechies 3)', avg_mse_W1, var(mse_W1_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W2 (Daubechies 3)', avg_mse_W2, var(mse_W2_vec));


 % Boxplot of Average MSEs
 figure;
 boxplot([mse_W1W2_vec, mse_W1_vec, mse_W2_vec], ...
     {'W1W2', 'Haar', 'Daubechies 3'});
 ylabel('MSE');
 title('MSE Distributions for M = 200 Simulations');

 figure;
 boxplot([mse_kron_vec, mse_W1W2W1_vec, mse_W1W2_vec, mse_W1_vec, mse_W2_vec], ...
     {'kron', 'W1W2W1', 'W1W2', 'Haar', 'Daubechies 3'});
 ylabel('MSE');
 title('MSE Distributions including Kron and W1W2W1');

%figure;
%histogram(mse_W1W2W1_vec, 15);
%title('Histogram of W1W2W1 MSEs');
%xlabel('MSE'); ylabel('Frequency');

%figure;
%histogram(mse_kron_vec, 15);
%title('Histogram of Kron MSEs');
