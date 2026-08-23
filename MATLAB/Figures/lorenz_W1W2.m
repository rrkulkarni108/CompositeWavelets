%% Lorenz curves: W1W2 vs W1 vs W2 on DJ signals
% Uses power ratio definition of SNR dB
close all; clear; clc;
rng(1);         

% User defined parameters
M      = 200;   % number of simulations
SNR_dB = 5;     % Power ratio in dB: 10*log10(P_signal/P_noise),
                % P_signal = ||s||^2/N, P_noise = sigma^2 = 1
n      = 1024;
signals = {'Doppler','Blocks','Bumps','HeaviSine'};  % order of subplots

% colors
cW1W2 = [0.00 0.45 0.70];   % blue
cW1   = [0.80 0.20 0.20];   % red
cW2   = [1.00, 0.50, 0.00];

% filters
haar = [1/sqrt(2), 1/sqrt(2)];
sym4 = [-0.075765714789273, -0.029635527645999, 0.497618667632015, ...
         0.803738751805916,  0.297857795605605, -0.099219543576935, ...
        -0.012603967262038,  0.032223100604043];
db3  = [0.332670552950165, 0.806891509311400, 0.459877502118228, ...
       -0.135011020010067, -0.085441273882042, 0.035226291882017];
db5  = [ 0.0033357252854737713, -0.012580751999015526, -0.0062414902130117052, ...
         0.077571493840045713, -0.032244869585029520, -0.24229488706619015, ...
         0.13842814590110342, 0.72430852843857444, 0.60382926979718952, ...
         0.16010239797419290];

% create a structure to hold signal computations
cfg = struct();

% Choose the signals which are most appropriate for the DJ signal (based on
% Donoho Johnstone 1995)

% Doppler: W1 = Daub5, W2 = Symm4
cfg.Doppler.W1 = db5;   cfg.Doppler.W2 = sym4;   cfg.Doppler.nl = 3;

% Blocks: W1 = Haar, W2 = Haar
cfg.Blocks.W1 = haar;   cfg.Blocks.W2 = haar;    cfg.Blocks.nl = 3;

% Bumps: W1 = Haar, W2 = Daub3  (and nl=1)
cfg.Bumps.W1  = haar;   cfg.Bumps.W2  = db3;     cfg.Bumps.nl  = 1;

% HeaviSine: W1 = Symm4, W2 = Symm4
cfg.HeaviSine.W1 = sym4; cfg.HeaviSine.W2 = sym4; cfg.HeaviSine.nl = 3;

% Helper function to calculate Lorenz from vector of energies
% given vector v, returns normalized Lorenz L(p) with length n (p in [0,1])
lorenz_of = @(v) ( ...
    cumsum(sort(abs(v(:)).^2,'ascend')) ./ sum(abs(v(:)).^2) );

% X-axis for plotting
p = linspace(0,1,n)';           % p in [0,1]
x_percent = p * 100;

%  Figure 2×2 subplots of Lorenz curves
fig = figure('Color','w','Units','inches','Position',[1 1 7.5 5.6], ...
    'Name','Lorenz curves: W1W2 vs W1 vs W2','NumberTitle','off');

for si = 1:numel(signals)
    name = signals{si};
    nl_use = cfg.(name).nl;

    % first build clean signal via MakeSignal then scale to the 
    % Power ratio (second moment mean(s.^2), instead of previous var) so that
    % P_signal/P_noise = 10^(SNR_dB/10) with unit-variance noise
    s_clean = MakeSignal(name, n);
    s_clean = s_clean(:);
    s_clean = s_clean ./ sqrt(mean(s_clean.^2)) * sqrt(10^(SNR_dB/10));

    % Wavelet matrices
    W1  = Wavmat(cfg.(name).W1, n, nl_use, 0);
    W2  = Wavmat(cfg.(name).W2, n, nl_use, 0);
    W12 = W1 * W2;

    includeW2 = ~isequal(cfg.(name).W1(:), cfg.(name).W2(:));

    % preallocate the values for lorenz into vectors
    L_W1W2 = zeros(n, M);
    L_W1   = zeros(n, M);
    if includeW2, L_W2 = zeros(n, M); end

    % M Simulations
    for m = 1:M
        s_noisy = s_clean + randn(n,1);   % unit-variance noise => P_noise = 1

        % get coefficients
        c12 = W12 * s_noisy;     % W1W2
        c1  = W1  * s_noisy;     % W1
        if includeW2
            c2 = W2 * s_noisy;   % W2
        end

        % Normalized Lorenz curves
        L_W1W2(:,m) = lorenz_of(c12);
        L_W1(:,m)   = lorenz_of(c1);
        if includeW2
            L_W2(:,m) = lorenz_of(c2);
        end
    end

    % Mean curves over the M simulations
    mu_W1W2 = mean(L_W1W2, 2);
    mu_W1   = mean(L_W1,   2);
    if includeW2
        mu_W2 = mean(L_W2, 2);
    end

    % Panel
    ax = subplot(2,2,si); hold on; grid on;

    plot(x_percent, mu_W1W2, 'Color', cW1W2, 'LineWidth', 1.3);
    plot(x_percent, mu_W1,   'Color', cW1,   'LineWidth', 1.3);
    if includeW2
        plot(x_percent, mu_W2, 'Color', cW2, 'LineWidth', 1.3);
    end

    title(sprintf('%s (nl = %d)', name, nl_use), 'FontSize', 16, 'FontWeight','bold');
    xlabel('p% of sorted coefficients', 'FontSize', 14);
    ylabel('Lorenz L(p)', 'FontSize', 14);
    xlim([-5, 105]); ylim([0, 1.02]);
    set(ax, 'FontSize', 13, 'LineWidth', 1.4, ...
            'XColor','k', 'YColor','k', 'GridAlpha', 0.15, 'Box','on');

    if includeW2
        legend({'W1W2','W1','W2'}, 'Location','northwest', ...
               'Box','off', 'FontSize', 13);
    else
        legend({'W1W2','W1'}, 'Location','northwest', ...
               'Box','off', 'FontSize', 13);
    end

    hold off;
end

set(fig,'PaperUnits','inches','PaperPosition',[0 0 7.5 5.6],'PaperSize',[7.5 5.6]);
print(fig,'Figures/lorenz_W1W2.pdf','-dpdf','-r300');
