%% DOPPLER signal: 200 sims, W1 = complex Daub6, W2 = Haar  
% Comparison of product transform, similarity, Kronecker and single base
% transforms when one of the bases is complex
% Needs: Wavmat.m

close all
clear; clc;

M    = 200;        % number of simulations
N    = 1024;
L    = 3;          % levels
n    = N;
SNR  = 5;          % can be 3, 5 or 7

t      = linspace(0,1,N).';
s_orig = sqrt(t.*(1-t)) .* sin((2*pi*1.05)./(t+.05));   % Doppler 

% W1: Complex DAUB6
hfilt_c = [ -0.066291260736239 - 0.085581649610182i, ...
             0.110485434560398 - 0.085581649610182i, ...
             0.662912607362388 + 0.171163299220364i, ...
             0.662912607362388 + 0.171163299220364i, ...
             0.110485434560398 - 0.085581649610182i, ...
            -0.066291260736239 - 0.085581649610182i ];

% W2: Haar
hfilt_r = [1/sqrt(2), 1/sqrt(2)];

W1   = Wavmat(hfilt_c, N, L);          % complex Daub6
W2   = Wavmat(hfilt_r, N, L);          % Haar
W12  = W1 * W2;                        % product 
W121 = W1' * W2 * W1;                  % similarity transform 

% Kronecker exploration 
%W0    = Wavmat([0.332670552950165  0.806891509311400  0.459877502118228 ...
              % -0.135011020010067 -0.085441273882042  0.035226291882017], 2^7, 4, 0); % Daub3, 128x128
%W0_2  = Wavmat([1/sqrt(2) 1/sqrt(2)], 2^3, 1, 0);                                      % Haar,   8x8
%WW    = kron(W0, W0_2);               % 1024 x 1024

%  Kronecker product with W1 = complex Daub6 and W2 = real Haar
% choose sizes whose product is N=1024; 128 x 8
N1 = 2^7;  
L1 = min(L, floor(log2(N1)));   % 128
N2 = 2^3;  
L2 = min(L, floor(log2(N2)));   % 8

Wc_small = Wavmat(hfilt_c, N1, L1);        % complex Daub6 
Wr_small = Wavmat(hfilt_r, N2, L2);        % real Haar 

WW = kron(Wc_small, Wr_small);             


% Universal hard threshold 
lambda = sqrt(2 * log(n));

mse_kron_vec   = zeros(M,1);
mse_W1W2W1_vec = zeros(M,1);
mse_W1W2_vec   = zeros(M,1);
mse_W1_vec     = zeros(M,1);
mse_W2_vec     = zeros(M,1);

%  Simulations
for m = 1:M
    %  SNR scaling
    s_base = s_orig;                                
    signal_std = sqrt(var(s_base));                  % std 
    s_snr = (s_base / signal_std) * sqrt(SNR);       % std(s_snr) = sqrt(SNR)
    s = s_snr + randn(N,1);                          % add N(0,1) noise

    %%  W1' W2 W1 (similarity)
    wd121 = W121 * s;
    wd121(abs(wd121) < lambda) = 0;
    sr121 = W121' * wd121;
    mse_W1W2W1_vec(m) = mean((s_snr - real(sr121)).^2);

    %%  W1 W2 product
    wd12 = W12 * s;
    wd12(abs(wd12) < lambda) = 0;
    sr12 = W12' * wd12;
    mse_W1W2_vec(m) = mean((s_snr - real(sr12)).^2);

    %%  Single base: W1 (complex Daub6)
    wd1 = W1 * s;
    wd1(abs(wd1) < lambda) = 0;
    sr1 = W1' * wd1;
    mse_W1_vec(m) = mean((s_snr - real(sr1)).^2);

    %%  Single base: W2 (Haar)
    wd2 = W2 * s;
    wd2(abs(wd2) < lambda) = 0;
    sr2 = W2' * wd2;
    mse_W2_vec(m) = mean((s_snr - real(sr2)).^2);

    %%  Kronecker product
    wd_kron = WW * s;
    wd_kron(abs(wd_kron) < lambda) = 0;
    sr_kron = WW' * wd_kron;
    mse_kron_vec(m) = mean((s_snr - real(sr_kron)).^2);
end

%  Summary 
avg_mse_kron   = mean(mse_kron_vec);
avg_mse_W1W2W1 = mean(mse_W1W2W1_vec);
avg_mse_W1W2   = mean(mse_W1W2_vec);
avg_mse_W1     = mean(mse_W1_vec);
avg_mse_W2     = mean(mse_W2_vec);

fprintf('| %-22s | %-11s | %-16s |\n', 'Method', 'Average MSE', 'Variance of MSE');
fprintf('|------------------------|-------------|------------------|\n');
fprintf('| %-22s | %.4f      | %.4f           |\n', 'Kronecker Product', avg_mse_kron,   var(mse_kron_vec));
fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1''W2W1 Product',  avg_mse_W1W2W1, var(mse_W1W2W1_vec));
fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1W2 Product',      avg_mse_W1W2,   var(mse_W1W2_vec));
fprintf('| %-22s | %.4f      | %.4f           |\n', 'W1 (complex D6)',   avg_mse_W1,     var(mse_W1_vec));
fprintf('| %-22s | %.4f      | %.4f           |\n', 'W2 (Haar)',         avg_mse_W2,     var(mse_W2_vec));

%  Boxplots 
figure;
boxplot([mse_W1W2_vec, mse_W1_vec, mse_W2_vec], ...
        {'W1W2','W1 (cD6)','W2 (Haar)'});
ylabel('MSE'); title(sprintf('Doppler: MSE over M = %d simulations (SNR=%g)', M, SNR));
colors = [0.8 0.3 0.3; 0.2 0.6 1.0; 0.2 0.6 1.0];
boxes = findobj(gca,'Tag','Box');
for j = 1:length(boxes)
    patch(get(boxes(j),'XData'), get(boxes(j),'YData'), colors(length(boxes)-j+1,:), ...
        'FaceAlpha',0.5, 'EdgeColor','k');
end
figure;
boxplot([mse_kron_vec, mse_W1W2W1_vec, mse_W1W2_vec, mse_W1_vec, mse_W2_vec], ...
        {'kron','W1''W2W1','W1W2','W1 (cD6)','W2 (Haar)'});
ylabel('MSE'); title(sprintf('Doppler: MSE Distributions , M = %d, SNR=%g', M, SNR));

colors = [0.6 0.2 0.9; 0.6 0.2 0.9; 0.6 0.2 0.9;  0.2 0.6 1.0; 0.2 0.6 1.0];
boxes = findobj(gca,'Tag','Box');
for j = 1:length(boxes)
    patch(get(boxes(j),'XData'), get(boxes(j),'YData'), colors(length(boxes)-j+1,:), ...
        'FaceAlpha',0.5, 'EdgeColor','k');
end

