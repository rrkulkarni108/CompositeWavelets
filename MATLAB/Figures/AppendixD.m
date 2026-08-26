%% Multiple SNR robustness on 4 DJ signals: product W1W2 vs single W1
%  In Appendix D
%  SNR is Power ratio in dB:
%      SNR = 10*log10( P_signal / P_noise ), P = mean(x.^2) aka P_signal = ||x||^2/N,
%      P_noise = sigma^2.
%
%  Convention is that sigma is fixed at 1 and the Signal is rescaled to hit each
%  target dB. This keeps the universal threshold at sigma*sqrt(2 log N)
%  with sigma = 1 and makes the 5 dB column directly comparable with the
%  main simulation results (Figure 1, Figure 2, etc).

close all; clear; clc; rng(0);

n      = 1024;
M      = 200;
nl     = 3;                        % Bumps overridden to 1 below
sigma  = 1;                        % noise sd, fixed
lambda = sigma * sqrt(2*log(n));   % universal threshold
SNRdB  = [-5 0 5 8.45];           % range of SNR dBs to empirically check different noise regimes

signals = {'Doppler','Blocks','Bumps','HeaviSine'};

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

cfg = struct();
cfg.Doppler.W1   = db5;  cfg.Doppler.W2   = sym4;
cfg.Blocks.W1    = haar; cfg.Blocks.W2    = haar;
cfg.Bumps.W1     = haar; cfg.Bumps.W2     = db3;
cfg.HeaviSine.W1 = sym4; cfg.HeaviSine.W2 = sym4;

mse_prod = zeros(numel(signals), numel(SNRdB));
mse_sing = zeros(numel(signals), numel(SNRdB));

for si = 1:numel(signals)

    name = signals{si};
    if strcmp(name,'Bumps'), L = 1; else, L = nl; end

    W1 = Wavmat(cfg.(name).W1, n, L, 0);
    W2 = Wavmat(cfg.(name).W2, n, L, 0);

    x0 = MakeSignal(name, n);
    x0 = x0(:);                                 

    % orthogonality check
    W12 = W1 * W2;
    if abs(norm(W12*x0)^2 - norm(x0)^2) > 1e-8 * norm(x0)^2
        warning('%s: energy not preserved -- check Wavmat', name);
    end

    for k = 1:numel(SNRdB)

        % rescale Signal to the target power ratio 
        x = x0 ./ sqrt(mean(x0.^2)) * sqrt(sigma^2 * 10^(SNRdB(k)/10));

        accP = 0; accS = 0;
        for mrep = 1:M
            y = x + sigma*randn(n,1);            % same for both transforms

            % product W1*W2
            cP = W12 * y;
            cP(abs(cP) < lambda) = 0;
            srP = W12' * cP;
            accP = accP + mean((x - srP).^2);

            % single W1
            cS = W1 * y;
            cS(abs(cS) < lambda) = 0;
            srS = W1' * cS;
            accS = accS + mean((x - srS).^2);
        end

        mse_prod(si,k) = accP / M;
        mse_sing(si,k) = accS / M;
    end
end

ratio = mse_sing ./ mse_prod;     % ratio of single to prod, if >1 then prod is better than single

%% Figure 10 -- AMSE vs SNR_in
fig = figure('Color','w','Units','inches','Position',[1 1 7.0 5.2]);
t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
t.OuterPosition = [0 0 1 0.93];

hLines = gobjects(1,2);
for si = 1:numel(signals)
    ax = nexttile; hold(ax,'on');
    p1 = plot(SNRdB, mse_prod(si,:), '-o','LineWidth',2.0,'MarkerSize',9,...
        'MarkerFaceColor','w');
    p2 = plot(SNRdB, mse_sing(si,:), '--s','LineWidth',2.0,'MarkerSize',9,...
        'MarkerFaceColor','w');
    if si == 1, hLines(1) = p1; hLines(2) = p2; end

    title(signals{si},'FontSize',20,'FontWeight','bold');
    xlabel('SNR_{in} (dB)','FontSize',16);
    ylabel('AMSE','FontSize',16);
    set(ax,'FontSize',15,'LineWidth',1.4,'XTick',SNRdB,...
        'XLim',[min(SNRdB)-1.5, max(SNRdB)+1.5],...
        'XColor','k','YColor','k','Box','on');
    yl = ylim(ax);
    ylim(ax, [yl(1) - 0.02*diff(yl), yl(2) + 0.06*diff(yl)]);
end

lg = legend(hLines, {'Product W_1W_2','Single W_1'}, ...
    'Orientation','vertical','FontSize',14,'Box','on',...
    'EdgeColor','k','Color','w');
lg.Units = 'normalized';
drawnow;
lg.Position = [0.5-lg.Position(3)/2, 0.44, lg.Position(3), lg.Position(4)];

axT = axes(fig,'Position',[0 0.94 1 0.05],'Visible','off');
text(axT, 0.5, 0.5, 'Denoising MSE vs Noise Level', ...
    'FontSize',17,'FontWeight','bold', ...
    'HorizontalAlignment','center','VerticalAlignment','middle');

drawnow;
set(fig,'PaperUnits','inches','PaperPosition',[0 0 7.0 5.2]);
print(fig,'Figures/mse_vs_snr.png','-dpng','-r300');
% output
fprintf('\nSNR_in (dB):   '); fprintf('%8d ', SNRdB); fprintf('\n');
for si = 1:numel(signals)
    fprintf('\n%s\n', signals{si});
    fprintf('  product      '); fprintf('%8.4f ', mse_prod(si,:)); fprintf('\n');
    fprintf('  single       '); fprintf('%8.4f ', mse_sing(si,:)); fprintf('\n');
    fprintf('  ratio S/P    '); fprintf('%8.2f ', ratio(si,:));    fprintf('\n');
end
