
%% Adaptive case - many (M) number of simulations - using SNR
% Comparison of single base wavelet matrix with adaptive method

close all
clear; clc;

% User Defined Parameters 
M = 200;  % number of simulations
SNR = 5;  % can change to 3, 5, or 7 


% Storing MSE values
mse_adaptive_vec = zeros(M, 1); 
mse_W1_vec = zeros(M, 1);
mse_W2_vec = zeros(M, 1);
mse_W3_vec = zeros(M, 1);
mse_W4_vec = zeros(M, 1);


n1=512;
n=1024 *2 ;

% Signal 1: Doppler
m=n1;
t = (1:m) ./m;
sig = sqrt(t.*(1-t)).*sin((2*pi*1.05) ./(t+.05));
s1=sig;


% Signal 2: Blocks
m=n1;
t = (1:m) ./m;
pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
hgt = [4 (-5) 3 (-4) 5 (-4.2) 2.1 4.3  (-3.1) 2.1 (-4.2)];
sig = zeros(size(t));
for j=1:length(pos)
    sig = sig + (1 + sign(t-pos(j))).*(hgt(j)/2) ;
end
s2=sig;


% Signal 3 : HeaviSine
m=n1;
t = (1:m) ./m;
sig = 4.*sin(4*pi.*t);
sig = sig - sign(t - .3) - sign(.72 - t);
s3 = sig;


% Signal 4: Bumps
m=n1;
t = (1:m) ./m;
pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
hgt = [ 4  5   3   4  5  4.2 2.1 4.3  3.1 5.1 4.2];
wth = [.005 .005 .006 .01 .01 .03 .01 .01  .005 .008 .005];
sig = zeros(size(t));
for j =1:length(pos)
   sig = sig + hgt(j)./( 1 + abs((t - pos(j))./wth(j))).^4;
end 
s4=sig;

s_orig = [s1 0.2 * s2 0.1 * s3 0.2 * s4]; % combined signal

I1=zeros(n);
for i=1:n1 
    I1(i,i)=1;
end


I2=zeros(n);
for i=n1+1:2*n1
    I2(i,i)=1;
end

 I3=zeros(n);
for i=2*n1+1:3*n1
    I3(i,i)=1;
end


 I4=zeros(n);
for i=3*n1+1:n
    I4(i,i)=1;
end

% Filter vectors 
hfilt1 = [ -0.07576571478934  -0.02963552764595  ...
          0.49761866763246   0.80373875180522  ...
          0.29785779560554  -0.09921954357694  ...
         -0.01260396726226   0.03222310060407];  % Symmlet 4 for Doppler

 hfilt2=[sqrt(2)/2 sqrt(2)/2]; % Haar for Blocks

 hfilt3=[0.230377813308959   0.714846570552874   0.630880767929889  -0.027983769416995 ...
     -0.187034811718967 0.030841381835995   0.032883011666994  -0.010597401784998]; % Daub4 for HeaviSine


 hfilt4 = [0.332670552950165   0.806891509311400   0.459877502118228  ...
     -0.135011020010067  -0.085441273882042  0.035226291882017]; % Daub3 for Bumps

 nl=3;        
 W1= Wavmat(hfilt1, n, nl, 4); % Wavelet matrix 
 W2= Wavmat(hfilt2, n, nl, 0); 
 W3= Wavmat(hfilt3, n, nl, 2);  
 W4= Wavmat(hfilt4, n, nl, 2);  


 W = W1 * I1 + W2 * I2 + W3 * I3 + W4 * I4;  % combine filters

 for m = 1:M 
     % Add Gaussian noise with fixed SNR (user can set at beginning of file)
     signal_var = var(s_orig);
     signal_std = sqrt(signal_var);
     s_orig= s_orig/signal_std * sqrt(SNR); % signal standardized to have unit variance, then scaled by sqrt(SNR)
     s = s_orig + randn(size(s_orig)); % adds unit variance Gaussian noise

     %%  Adaptive Wavelet Shrinkage 
     wd = W * s'; % apply wavelet to noisy signal

     n = length(s);
     lambda = sqrt(2 * log(n)); % Universal thresholding, noise ~ Normal(0,1)
    
     wd(abs(wd) < lambda) = 0;  % Hard thresholding
    
     sr = W' * wd; % reconstruct denoised signal

     mse_adaptive_vec(m) = mean((s_orig - sr').^2); 


     %%  Single Base (W1 - Symmlet 4) Wavelet Shrinkage  

     wd_W1 = W1 * s'; % Apply single wavelet matrix W1 to noisy signal

     wd_thresh_W1 = wd_W1;
     wd_thresh_W1(abs(wd_W1) < lambda) = 0; % Hard thresolding
    
     sr_W1 = W1' * wd_thresh_W1; % reconstruct denoised signal
    
     mse_W1_vec(m) = mean((s_orig' - sr_W1).^2); % mse for single wavelet W1
    

     %% Single Base (W2 - Haar) Wavelet Shrinkage  

     wd_W2 = W2 * s'; % Apply single wavelet matrix W2 to noisy signal

     wd_thresh_W2 = wd_W2;
     wd_thresh_W2(abs(wd_W2) < lambda) = 0; % Hard thresolding
    
     sr_W2 = W2' * wd_thresh_W2; % reconstruct denoised signal
    
     mse_W2_vec(m) = mean((s_orig' - sr_W2).^2); % mse for single wavelet W2
     
     
     %% Single Base (W3 - DAUB4) Wavelet Shrinkage  

     wd_W3 = W3 * s'; % Apply single wavelet matrix W3 to noisy signal

     wd_thresh_W3 = wd_W3;
     wd_thresh_W3(abs(wd_W3) < lambda) = 0; % Hard thresolding
    
     sr_W3 = W3' * wd_thresh_W3; % reconstruct denoised signal
    
     mse_W3_vec(m) = mean((s_orig' - sr_W3).^2); % mse for single wavelet W3
    
     %% Single Base (W4 - DAUB 3) Wavelet Shrinkage  

     wd_W4 = W4 * s'; % Apply single wavelet matrix W4 to noisy signal

     wd_thresh_W4 = wd_W4;
     wd_thresh_W4(abs(wd_W4) < lambda) = 0; % Hard thresolding
    
     sr_W4 = W4' * wd_thresh_W4; % reconstruct denoised signal
    
     mse_W4_vec(m) = mean((s_orig' - sr_W4).^2); % mse for single wavelet W4
    

 end

 % Calculate Average MSE across M simulations for Adaptive W
 % and all 4 single bases W1-W4
 avg_mse_adaptive = mean(mse_adaptive_vec);
 avg_mse_W1 = mean(mse_W1_vec);
 avg_mse_W2 = mean(mse_W2_vec);
 avg_mse_W3 = mean(mse_W3_vec);
 avg_mse_W4 = mean(mse_W4_vec);

 % Print results
 fprintf('Average MSE (Adaptive): %.4f\n', avg_mse_adaptive);
 fprintf('Average MSE (W1 - Sym4):   %.4f\n', avg_mse_W1);
 fprintf('Average MSE (W2 - Haar):   %.4f\n', avg_mse_W2);
 fprintf('Average MSE (W3 - Daub4):   %.4f\n', avg_mse_W3);
 fprintf('Average MSE (W4 - Daub3):   %.4f\n', avg_mse_W4);

 % variances
 var_mse_adaptive = var(mse_adaptive_vec);
 var_mse_W1 = var(mse_W1_vec);
 var_mse_W2 = var(mse_W2_vec);
 var_mse_W3 = var(mse_W3_vec);
 var_mse_W4 = var(mse_W4_vec);

 % Boxplot of Average MSEs
 figure;
 boxplot([mse_adaptive_vec, mse_W1_vec, mse_W2_vec, mse_W3_vec, mse_W4_vec], ...
     {'Adaptive', 'Symmlet 4', 'Haar', 'DAUB 4', 'DAUB 3'});
 ylabel('MSE');
 title('MSE Distributions');


colors = [0.6 0.2 0.9;  0.2 0.6 1.0;  0.2 0.6 1.0; 0.2 0.6 1.0; 0.2 0.6 1.0];
boxes = findobj(gca,'Tag','Box');
for j = 1:length(boxes)
    patch(get(boxes(j),'XData'), get(boxes(j),'YData'), colors(length(boxes)-j+1,:), ...
        'FaceAlpha',0.5, 'EdgeColor','k');
end

% Printing table in command window
fprintf('| %-12s | %-10s | %-10s |\n', 'Method', 'Average MSE', 'Variance');
fprintf('|-------------|------------|------------|\n');
fprintf('| %-12s | %.4f     | %.4f     |\n', 'Adaptive', avg_mse_adaptive, var_mse_adaptive);
fprintf('| %-12s | %.4f     | %.4f     |\n', 'Symmlet 4', avg_mse_W1, var_mse_W1);
fprintf('| %-12s | %.4f     | %.4f     |\n', 'Haar', avg_mse_W2, var_mse_W2);
fprintf('| %-12s | %.4f     | %.4f     |\n', 'Daub 4', avg_mse_W3, var_mse_W3);
fprintf('| %-12s | %.4f     | %.4f     |\n', 'Daub 3', avg_mse_W4, var_mse_W4);
