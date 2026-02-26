%% Lorenz curves: W1W2 vs W1 vs W2 on DJ signals
close all; clear; clc;


% User defined parameters
M   = 200;      % number of simulations
SNR = 5;        %  3, 5, 7
n   = 1024;
signals = {'Doppler','Blocks','Bumps','HeaviSine'};  % order of subplots

% colors
cW1W2 = [0.00 0.45 0.70];   % blue 
cW1   = [0.80 0.20 0.20];   % red
cW2   = [1.00, 0.50, 0.00] %[0.00 0.62 0.45];   % purple
bandAlpha = 0.15;
lineW     = 1;

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
    cumsum(sort(v(:).^2,'ascend')) ./ sum(v(:).^2) );

% X-axis for plotting
p = linspace(0,1,n)';           % p in [0,1]
x_percent = p * 100;

%  Figure 2×2 subplots of Lorenz curves
figure('Color','w','Name','Lorenz curves: W1W2 vs W1 vs W2','NumberTitle','off');

% To harmonize y-axis across panels if you like
globalY = [];

for si = 1:numel(signals)
    name = signals{si};
    nl_use = cfg.(name).nl;

    % first build clean signal via MakeSignal then scale using SNR
    s_clean = MakeSignal(name, n);
    s_clean = s_clean ./ std(s_clean) * sqrt(SNR);

    % Wavelet matrices  
    W1 = Wavmat(cfg.(name).W1, n, nl_use, 0);
    W2 = Wavmat(cfg.(name).W2, n, nl_use, 0);
    W12 = W1 * W2;   

    
    includeW2 = ~isequal(cfg.(name).W1(:), cfg.(name).W2(:));

    % preallocate the values for lorenz into vectors
    L_W1W2 = zeros(n, M);
    L_W1   = zeros(n, M);
    if includeW2, L_W2 = zeros(n, M); end

    % M Simulations
    for m = 1:M
        s_noisy = s_clean + randn(size(s_clean));  % unit-variance noise

        % get coefficents
        c12 = W12 * s_noisy';     % W1W2
        c1  = W1  * s_noisy';     % W1
        if includeW2
            c2 = W2  * s_noisy';  % W2
        end

        % Normalized Lorenz curves
        L_W1W2(:,m) = lorenz_of(c12);
        L_W1(:,m)   = lorenz_of(c1);
        if includeW2
            L_W2(:,m) = lorenz_of(c2);
        end
    end

    % Mean and +/-1 SD intervals
    mu_W1W2 = mean(L_W1W2, 2);  s
    d_W1W2 = std(L_W1W2, 0, 2);
    mu_W1   = mean(L_W1,   2);  s
    d_W1   = std(L_W1,   0, 2);
    if includeW2
        mu_W2 = mean(L_W2, 2);  sd_W2   = std(L_W2,   0, 2);
    end

    % Panel
    subplot(2,2,si); hold on; grid on;

    % Bands 
    % W1W2 band
    fill([x_percent; flipud(x_percent)], ...
         [mu_W1W2 - sd_W1W2; flipud(mu_W1W2 + sd_W1W2)], ...
         cW1W2, 'FaceAlpha', bandAlpha, 'EdgeColor', 'none');

    % W1 band
    fill([x_percent; flipud(x_percent)], ...
         [mu_W1 - sd_W1; flipud(mu_W1 + sd_W1)], ...
         cW1, 'FaceAlpha', bandAlpha, 'EdgeColor', 'none');

    % W2 band (if W2 != W1)
    if includeW2
        fill([x_percent; flipud(x_percent)], ...
             [mu_W2 - sd_W2; flipud(mu_W2 + sd_W2)], ...
             cW2, 'FaceAlpha', bandAlpha, 'EdgeColor', 'none');
    end

    % Mean curves
    plot(x_percent, mu_W1W2, 'Color', cW1W2, 'LineWidth', lineW);
    plot(x_percent, mu_W1,   'Color', cW1,   'LineWidth', lineW);
    if includeW2
        plot(x_percent, mu_W2, 'Color', cW2, 'LineWidth', lineW);
    end

    title(sprintf('%s (nl = %d)', name, nl_use));
    xlabel('Percent'); ylabel('Lorenz L(p)'); xlim([-5, 105]); ylim([0, 1.02]);

    globalY = [globalY; mu_W1W2; mu_W1];
    if includeW2, globalY = [globalY; mu_W2]; end

    % Legend per panel 
    if includeW2
        legend({'W1W2 band','W1 band','W2 band','W1W2','W1','W2'}, ...
               'Location','southeast','Box','off');
    else
        legend({'W1W2 band','W1 band','W1W2','W1'}, ...
               'Location','southeast','Box','off');
    end

    hold off;
end
