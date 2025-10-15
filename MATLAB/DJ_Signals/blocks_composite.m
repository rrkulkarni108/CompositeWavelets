%% BLOCKS from Donoho-Johnstone Library of Signals

%% Matrix Product case - many (M) number of simulations - using SNR
% Comparison of single base wavelet matrix with Matrix Products 
%% First: Two-Matrix Product W1W2 comparison to Single base
%% Second: Similarity Transform W1'W2W1 comparison to W1W2, Single base
%% Third: Kronecker Product 


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


% Signal: Doppler
%n = 1024;
%t = linspace(0, 1, n);
%doppler = sqrt(t .* (1 - t)) .* sin((2 * pi * 1.05) ./ (t + 0.05));
%s_orig = doppler;  

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


     %%  Single Base (W1 - Haar) Wavelet Shrinkage  

     wd_W1 = W1 * s'; % Apply single wavelet matrix W1 to noisy signal

     wd_thresh_W1 = wd_W1;
     wd_thresh_W1(abs(wd_W1) < lambda) = 0; % Hard thresolding
    
     sr_W1 = W1' * wd_thresh_W1; % reconstruct denoised signal
    
     mse_W1_vec(m) = mean((s_orig' - sr_W1).^2); % mse for single wavelet W1
    

     %%  Single Base (W2 - Haar) Wavelet Shrinkage  

     wd_W2 = W2 * s'; % Apply single wavelet matrix W2 to noisy signal

     wd_thresh_W2 = wd_W2;
     wd_thresh_W2(abs(wd_W2) < lambda) = 0; % Hard thresolding
    
     sr_W2 = W2' * wd_thresh_W2; % reconstruct denoised signal
    
     mse_W2_vec(m) = mean((s_orig' - sr_W2).^2); % mse for single wavelet W2
     
    %%  Kronecker Product 
    % Exploration of kron(W0,W0) 
     W0 = Wavmat(hfilt_haar, 2^7, 4, 0); % Wavelet matrix 
     %W0_2 = Wavmat(hfilt1, 2^3, 1, 0); % with level 1, variability is smaller
     W0_2 = Wavmat([sqrt(2)/2 sqrt(2)/2], 2^3, 1, 0); % with level 1, variability is smaller

     WW = kron(W0,W0_2);
     wd_kron = WW * s'; % apply wavelet to noisy signal
     wd_kron(abs(wd_kron) < lambda) = 0;  % Hard thresholding
        
     sr_kron = WW' * wd_kron; % reconstruct denoised signal
    
     mse_kron_vec(m) = mean((s_orig' - sr_kron).^2); % mse for kronecker prod


 end

 % Calculate Average MSE for kron(W0, W0_2) across M simulations 
 avg_mse_kron = mean(mse_kron_vec);

 % Calculate Average MSE for W1'W2W1 across M simulations 
 avg_mse_W1W2W1 = mean(mse_W1W2W1_vec);

 % Calculate Average MSE for W1W2 across M simulations 
 avg_mse_W1W2 = mean(mse_W1W2_vec);


 % and all 4 single bases W1-W4
 avg_mse_W1 = mean(mse_W1_vec);
 avg_mse_W2 = mean(mse_W2_vec);

 % Print results

 %fprintf( 'Average MSE (Kronecker Product):   %.4f\n', avg_mse_kron); 

 %fprintf( 'Variance of MSEs (Kronecker Product):   %.4f\n', var(mse_kron_vec)); 

 %fprintf( '-------------------------------------------\n'); 

 %fprintf('Average MSE (W1W2W1 Product): %.4f\n', avg_mse_W1W2W1);

 %fprintf('Average MSE (W1W2 Product): %.4f\n', avg_mse_W1W2);

 %fprintf('Average MSE (W1 - Haar):   %.4f\n', avg_mse_W1);
 %fprintf( 'Variance of MSEs (W1 - Haar):   %.4f\n', var(mse_W1_vec)); 
 %fprintf('Average MSE (W2 - Haar):   %.4f\n', avg_mse_W2);
 %fprintf( 'Variance of MSEs (W1 - Haar):   %.4f\n', var(mse_W1_vec)); 

 fprintf('| %-22s | %-11s | %-16s |\n', 'Method', 'Average MSE', 'Variance of MSE');
 fprintf('|------------------------|-------------|------------------|\n');
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'Kronecker Product', avg_mse_kron, var(mse_kron_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2W1 Product', avg_mse_W1W2W1, var(mse_W1W2W1_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2 Product', avg_mse_W1W2, var(mse_W1W2_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1 (Haar)', avg_mse_W1, var(mse_W1_vec));
 fprintf('| %-22s | %.4f      | %.4f           |\n', 'W2 (Haar)', avg_mse_W2, var(mse_W2_vec));


 % Boxplot of Average MSEs
 figure;
 boxplot([mse_W1W2_vec, mse_W1_vec, mse_W2_vec], ...
     {'W1W2', 'Haar', 'Haar'});
 ylabel('MSE');
 title('MSE Distributions for M = 200 Simulations');

 figure;
 boxplot([mse_kron_vec, mse_W1W2W1_vec, mse_W1W2_vec, mse_W1_vec, mse_W2_vec], ...
     {'kron', 'W1W2W1', 'W1W2', 'Haar', 'Haar'});
 ylabel('MSE');
 title('MSE Distributions including Kron and W1W2W1');

%figure;
%histogram(mse_W1W2W1_vec, 15);
%title('Histogram of W1W2W1 MSEs');
%xlabel('MSE'); ylabel('Frequency');

%figure;
%histogram(mse_kron_vec, 15);
%title('Histogram of Kron MSEs');
