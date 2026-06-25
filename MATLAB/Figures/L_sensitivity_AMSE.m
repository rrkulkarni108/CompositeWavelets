%% L-sensitivity: compressibility and MSE vs decomposition level
% For the four Donoho-Johnstone signals, composite W1W2 vs single basis.


close all; clear; clc; rng(0);

%  parameters 
n     = 1024;
SNR   = 5;            % variance ratio - var(x)/var(epsilon), SNR 5 = 7 dB.
M     = 200;          % Monte Carlo reps for MSE
Lvals = 1:5;          % decomposition levels L = 1, 2, 3, 4, 5.
signals = {'Doppler','Blocks','Bumps','HeaviSine'};

%  filters 
haar = [1/sqrt(2), 1/sqrt(2)];
sym4 = [-0.075765714789273,-0.029635527645999,0.497618667632015,...
    0.803738751805916,0.297857795605605,-0.099219543576935,...
    -0.012603967262038,0.032223100604043];
db3  = [0.332670552950165,0.806891509311400,0.459877502118228,...
    -0.135011020010067,-0.085441273882042,0.035226291882017];
db5  = [0.0033357252854737713,-0.012580751999015526,-0.0062414902130117052,...
    0.077571493840045713,-0.032244869585029520,-0.24229488706619015,...
    0.13842814590110342,0.72430852843857444,0.60382926979718952,...
    0.16010239797419290];

% W1 and W2 for each signal
cfg = struct();
cfg.Doppler.W1   = db5;  cfg.Doppler.W2   = sym4;
cfg.Blocks.W1    = haar; cfg.Blocks.W2    = haar;
cfg.Bumps.W1     = haar; cfg.Bumps.W2     = db3;
cfg.HeaviSine.W1 = sym4; cfg.HeaviSine.W2 = sym4;

% function that calculates top-k% energy concentration
topfrac = 0.05;   % top 5% of coefficients
topk_energy = @(c) local_topk(c, topfrac);
function f = local_topk(c, topfrac)
    e = sort(c(:).^2,'descend');
    k = ceil(topfrac*numel(e));
    f = sum(e(1:k))/sum(c(:).^2);
end

% storage of compression and mse in vectors
comp = zeros(numel(signals), numel(Lvals));   % compressibility of clean signal
mse  = zeros(numel(signals), numel(Lvals));   % AMSE after denoising


for si = 1:numel(signals)
    name = signals{si};
    % match Figs 1/2 in our paper: build using MakeSignal, SNR-scale 
    s0 = MakeSignal(name, n);
    s0 = s0 ./ std(s0) * sqrt(SNR);      % row vector

    for li = 1:numel(Lvals)
        L = Lvals(li);
        W1 = Wavmat(cfg.(name).W1, n, L, 0);
        W2 = Wavmat(cfg.(name).W2, n, L, 0);

        %  disbalance of clean signal
        c_clean = W1 * (W2 * s0');        % apply W2 then W1, on column s0'
        comp(si,li) = topk_energy(c_clean);

        %  MSE,  threshold + product + reconstruct
        lambda = sqrt(2*log(n));          % universal threshold, N(0,1) noise
        acc = 0;
        for m = 1:M
            s_noisy = s0 + randn(size(s0));      % row, unit-variance noise
            c = W1 * (W2 * s_noisy');            % product transform
            c(abs(c) < lambda) = 0;              % hard threshold 
            sr = W2' * (W1' * c);                % reconstruct
            acc = acc + mean((s0 - sr').^2);     % AMSE
        end
        mse(si,li) = acc / M;
    end
end



% First plot the compressibility vs L 
figure('Color','w'); hold on; box on;
markers = {'-o','-s','-d','-^'};
for si = 1:numel(signals)
    plot(Lvals, comp(si,:), markers{si}, 'LineWidth',1.6,'MarkerSize',7,...
        'DisplayName', signals{si});
end
xlabel('Decomposition level $L$','Interpreter','latex','FontSize',14);
ylabel('Top-5\% energy fraction','Interpreter','latex','FontSize',14);
set(gca,'FontSize',12,'XTick',Lvals);
legend('Location','southeast','FontSize',11);
title('Compressibility of $W_1 W_2$ vs level','Interpreter','latex','FontSize',14);

% Second plot the MSE vs L 
figure('Color','w'); hold on; box on;
for si = 1:numel(signals)
    plot(Lvals, mse(si,:), markers{si}, 'LineWidth',1.6,'MarkerSize',7,...
        'DisplayName', signals{si});
end
xlabel('Decomposition level $L$','Interpreter','latex','FontSize',14);
ylabel('AMSE after denoising','Interpreter','latex','FontSize',14);
set(gca,'FontSize',12,'XTick',Lvals);
legend('Location','northwest','FontSize',11);
title('Denoising MSE vs level','Interpreter','latex','FontSize',14);

% Table for printing
disp('Compressibility (top-5% energy) by signal x L:');
disp(array2table(comp,'VariableNames',compose('L%d',Lvals),'RowNames',signals));
disp('AMSE by signal x L:');
disp(array2table(mse, 'VariableNames',compose('L%d',Lvals),'RowNames',signals));
