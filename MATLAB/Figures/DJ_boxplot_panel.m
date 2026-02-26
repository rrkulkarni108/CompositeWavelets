

%% DJ 4-signal, 2x2 panels  -  W1W2, Kron, W1, (W2 if W1 != W2)
close all; clear; clc;

% User defined parameters
M   = 200;           % number of simulations
SNR = 5;             % can choose 3, 5, or 7
n   = 1024;
nl  = 3;             % WT levels 
lambda = sqrt(2*log(n));   % universal threshold for N(0,1) noise

signals = {'Doppler','Blocks','Bumps','HeaviSine'};  % order of the subplots

% Colors in order: W1W2, Kron, W1, W2
palette = [0 0.45 0.70; 0.835 0.369 0; 0 0.62 0.451; 0.8 0.475 0.655];

% Helper: build signal via MakeSignal and SNR-scale
mk_sig = @(name) (MakeSignal(name,n) ./ std(MakeSignal(name,n))) * sqrt(SNR);

% wavelet filters
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

% Map per signal: W1, W2, and Kronecker (W0 main, W0_2 small Haar)
cfg = struct();

% Doppler: W1=Daub5, W2=Symm4; Kron uses W0=db5(2^7,4), W0_2=Haar(2^3,1)
cfg.Doppler.W1  = db5;    cfg.Doppler.W2  = sym4;
cfg.Doppler.W0  = db5;    cfg.Doppler.W0_2 = haar;

% Blocks: W1=Haar, W2=Haar; Kron uses W0=Haar(2^7,4), W0_2=Haar(2^3,1)
cfg.Blocks.W1   = haar;   cfg.Blocks.W2   = haar;
cfg.Blocks.W0   = haar;   cfg.Blocks.W0_2 = haar;

% Bumps: W1=Haar, W2=Daub3; Kron uses W0=db3(2^7,4), W0_2=Haar(2^3,1)
cfg.Bumps.W1    = haar;   cfg.Bumps.W2    = db3;
cfg.Bumps.W0    = db3;    cfg.Bumps.W0_2  = haar;

% HeaviSine: W1=Symm4, W2=Symm4; Kron uses W0=sym4(2^7,4), W0_2=Haar(2^3,1)
cfg.HeaviSine.W1 = sym4;  cfg.HeaviSine.W2 = sym4;
cfg.HeaviSine.W0 = sym4;  cfg.HeaviSine.W0_2 = haar;

% Storage for global y-limits
allY = [];

% Make the boxplots figure with 2x2 subplots 
figure('Color','w'); 

for p = 1:4
    name = signals{p};
    s0   = mk_sig(name);  
    
    % choose nl per signal (Bumps gets 1, others keep global nl=3)
    if strcmp(name,'Bumps')
        nl_use = 1;
    else
        nl_use = nl;
    end 
    
    % Wavelet matrices for this signal
    W1 = Wavmat(cfg.(name).W1, n, nl_use , 0);
    W2 = Wavmat(cfg.(name).W2, n, nl_use , 0);
    
    % Kronecker basis for this signal
    W0   = Wavmat(cfg.(name).W0,   2^7, 4, 0);
    W0_2 = Wavmat(cfg.(name).W0_2, 2^3, 1, 0);
    WW   = kron(W0, W0_2);
    
    % omit W2 if W1==W2, otherwise include - flag
    includeW2 = ~isequal(cfg.(name).W1(:), cfg.(name).W2(:));
    
    % prespecify vectors for AMSE across M sims
    amse_W1W2 = zeros(M,1);
    amse_Kron = zeros(M,1);
    amse_W1   = zeros(M,1);
    if includeW2, amse_W2 = zeros(M,1); end
    
    % Simulation loop
    for m = 1:M
        s_noisy = s0 + randn(size(s0));   % unit-variance Gaussian noise
        
        % W1W2 product shrinkage
        c = (W1 * (W2 * s_noisy')); 
        c(abs(c) < lambda) = 0; % Hard thresholding
        sr = W2' * (W1' * c); % reconstruct denoised signal
        amse_W1W2(m) = mean((s0 - sr').^2); % calculate AMSE
        
        % Kronecker shrinkage
        ck = WW * s_noisy';
        ck(abs(ck) < lambda) = 0;
        srk = WW' * ck;
        amse_Kron(m) = mean((s0 - srk').^2);
        
        % W1 single-base shrinkage
        c1 = W1 * s_noisy'; 
        c1(abs(c1) < lambda) = 0;
        sr1 = W1' * c1;
        amse_W1(m) = mean((s0 - sr1').^2);
        
        % W2 single-base shrinkage (only if W1 != W2)
        if includeW2
            c2 = W2 * s_noisy'; 
            c2(abs(c2) < lambda) = 0;
            sr2 = W2' * c2;
            amse_W2(m) = mean((s0 - sr2').^2);
        end
    end
    
    % collect everything needed for global y-limits
    allY = [allY; amse_W1W2; amse_Kron; amse_W1];
    if includeW2, allY = [allY; amse_W2]; end
    
    %  here we begin to plot the whole panel
    subplot(2,2,p);
    hold on;
    
    % collect the data & group labels in desired order: W1W2, Kron, W1, (W2)
    if includeW2
        Y = [amse_W1W2; amse_Kron; amse_W1; amse_W2];
        G = [repmat("W1W2",M,1); repmat("Kron",M,1); repmat("W1",M,1); repmat("W2",M,1)];
        methodsPresent = {'W1W2','Kron','W1','W2'};
    else
        Y = [amse_W1W2; amse_Kron; amse_W1];
        G = [repmat("W1W2",M,1); repmat("Kron",M,1); repmat("W1",M,1)];
        methodsPresent = {'W1W2','Kron','W1'};
    end
    
    % Keep order
    G = categorical(G, methodsPresent, methodsPresent);
    
    % draw boxplot outline and fill with patches 
    boxplot(Y, G, 'Colors', repmat([0 0 0], numel(methodsPresent), 1), ...
                 'Symbol','+');
    title(name);
    ylabel('Average MSE');
    set(gca,'FontName','Times New Roman','FontSize',11);
    
    % fill boxes with colors
    bx = findobj(gca,'Tag','Box');
    bx = flipud(bx);
    for j = 1:numel(bx)
        set(bx(j), 'LineWidth', 0.8);
        patch(get(bx(j),'XData'), get(bx(j),'YData'), palette(j,:), ...
              'FaceAlpha',0.35, 'EdgeColor','k', 'LineWidth',0.8);
    end
    set(findobj(gca,'Tag','Median'),'Color',[0.15 0.15 0.15],'LineWidth',1.2);
    hold off;
end

%  align the y-lims across panels 
yl = [min(allY), max(allY)];
for p = 1:4
    subplot(2,2,p); ylim(yl);
end

% make one legend outside the grid 
presentUnion = {'W1W2','Kron','W1','W2'};
% only keep items that appear in at least one panel
keep = false(size(presentUnion));
for j = 1:numel(presentUnion)
    meth = presentUnion{j};
    for p = 1:4
        ax = subplot(2,2,p);
        cats = categories(ax.Children(1).XData); 
        if any(strcmp(cats, meth)), keep(j) = true; break; end
    end
end
presentUnion = presentUnion(keep);
ph = gobjects(numel(presentUnion),1);
figure; axL = axes('Visible','off'); hold(axL,'on');
for j = 1:numel(presentUnion)
    idx = find(strcmp({'W1W2','Kron','W1','W2'}, presentUnion{j}));
    ph(j) = plot(axL, nan, nan, 's', 'MarkerSize', 9, ...
                 'MarkerFaceColor', palette(idx,:), ...
                 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
end
legend(axL, ph, presentUnion, 'Orientation','horizontal','Box','off');
title(axL, 'Methods Legend');
