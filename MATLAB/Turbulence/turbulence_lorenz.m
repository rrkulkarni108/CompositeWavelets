% Turbulence data
% Lorenz curves for visual check of energy disbalance
%
% Columns in turbulence.txt:
%   1: U (streamwise velocity)
%   2: V (lateral velocity)
%   3: W (vertical velocity)
%   4: T (temperature in Kelvin)
%
% Wavelet bases:
%   W1 = Symmlet 4
%   W2 = Coiflet 3
%   Products: W1W2, W2W1
%
% We will analyze the first N = 2048 points of each component 
% as done in Katul and Vidakovic (1996).


clc; clear; close all;

data = load('turbulence.txt');

colNames = {'U (streamwise)', 'V (lateral)', ...
            'W (vertical)',   'T (temperature)'};

N = 2048;
if size(data,1) < N
    error('Not enough samples: need at least %d rows.', N);
end

% filters, same as Barbara
hfilt_symm4 = [-0.075765714789273, ...
               -0.029635527645999, ...
                0.497618667632015, ...
                0.803738751805916, ...
                0.297857795605605, ...
               -0.099219543576935, ...
               -0.012603967262038, ...
                0.032223100604043];

hfilt_coif3  = [ -0.003793512864,  0.007782596426,  0.023452696142, ...
                 -0.065771911281, -0.061123390003,  0.405176902410, ...
                  0.793777222626,  0.428483476378, -0.071799821619, ...
                 -0.082301927106,  0.034555027573,  0.015880544864, ...
                 -0.009007976137, -0.002574517688,  0.001117518771, ...
                  0.000466216960, -0.000070983303, -0.000034599773 ];

% Single and product bases
% Using the same last two parameters as in Barbara - 2,2
W1 = Wavmat(hfilt_symm4, N, 2, 2);   % Symm4
W2 = Wavmat(hfilt_coif3,  N, 2, 2);   % Coif3

W1W2 = W1 * W2;
W2W1 = W2 * W1;

transNames = {'Time domain', 'W1 (Symm4)', 'W2 (Coif3)', ...
              'W1W2',         'W2W1'};

%loop over U, V, W, T
for col = 1:4
    x = data(1:N, col);
    x = x(:);

    % Coefs in each domain
    c_time = x;
    c1     = W1   * x;
    c2     = W2   * x;
    c12    = W1W2 * x;
    c21    = W2W1 * x;

    coeffs = {c_time, c1, c2, c12, c21};

    % Energy summary (total, top 1%, top 5%) 
    fprintf('\n=============================================================\n');
    fprintf('Component: %s (column %d), N = %d samples\n', colNames{col}, col, N);
    fprintf('-------------------------------------------------------------\n');
    fprintf('Transform       TotalEnergy      Top1%%%%Energy   Top5%%%%Energy\n');
    fprintf('-------------------------------------------------------------\n');

    for t = 1:numel(coeffs)
        [E_total, frac1, frac5] = energy_summary(coeffs{t});
        fprintf('%-14s  %12.4f      %10.4f      %10.4f\n', ...
                transNames{t}, E_total, frac1, frac5);
    end

    % Lorenz curves for W only
    % col==1 if we want U
    if col == 3   % W (vertical)
        n = N;
        p = linspace(0,1,n);  % continuous parameter p in [0,1]
        floor_np = min(max(floor(n * p), 1), n); % floor(n p) clipped to [1,n]

        [lor_time] = lorenz_continuous(c_time, floor_np);

        % W1 Symm4
        [lor_W1]   = lorenz_continuous(c1, floor_np);

        % W2 Coif3
        [lor_W2]   = lorenz_continuous(c2, floor_np);

        % W1W2
        [lor_W1W2] = lorenz_continuous(c12, floor_np);

        % W2W1
        [lor_W2W1] = lorenz_continuous(c21, floor_np);

        %  Plot normalized Lorenz curves (L(p) in [0,1]) 
        figure;
        %plot(p, lor_time, 'k-',  'LineWidth', 1.5); 
        %plot(p, lor_W1,   'b-', 'LineWidth', 1.5); hold on;
        %plot(p, lor_W2,  '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
       % plot(p, lor_W1W2, 'r-',  'LineWidth', 1.5);
        %plot(p, lor_W2W1, 'k--',  'LineWidth', 1.5);
        % Single bases -lighter
        plot(p, lor_W1,  '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);hold on;
        plot(p, lor_W2,  '--','Color', [0 0.4470 0.7410], 'LineWidth', 2.5);
        
        % Product bases -bolder
        plot(p, lor_W1W2,'-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
        plot(p, lor_W2W1,'--','Color', [0.8500 0.3250 0.0980], 'LineWidth', 2.5);

        plot([0 1], [0 1], 'k:', 'LineWidth', 1); % equality line

        xlabel('p (proportion of coefficients)');
        ylabel('L(p) (cumulative energy)');
        title('Lorenz curves for W component (vertical velocity)');
        legend({ ...
            %'Time', ...
            'W1 (Symm4)', 'W2 (Coif3)', 'W1W2', 'W2W1', 'L(p)=p'}, ...
               'Location', 'northwest');
        xlim([-0.02, 1.02]); ylim([-0.02, 1.02]);
        grid on;
    end
end

% energy_summary (total energy, top 1%%, top 5%%)
function [E_total, frac1, frac5] = energy_summary(c)
    E_total = sum(c.^2);

    stats1 = energyFractions(c, 0.01);
    stats5 = energyFractions(c, 0.05);

    frac1 = stats1.E_top / E_total;
    frac5 = stats5.E_top / E_total;
end

% energyFractions (for a given proportion p)

function stats = energyFractions(c, p)
    E_total = sum(c.^2);
    k = max(1, round(p * numel(c)));
    e_sorted = sort(c.^2, 'descend');
    E_top = sum(e_sorted(1:k));

    stats.E_total = E_total;
    stats.E_top   = E_top;
end


%Lorenz curve (normalized)

function lorenz = lorenz_continuous(c, floor_np)
    % c: coefficient vector
    % floor_np: indices floor(n*p) for p in [0,1]
    e = c.^2;
    e_sorted = sort(e, 'ascend');    % Lorenz: sort from smallest to largest
    sums_e   = cumsum(e_sorted);     % cumulative energy
    totalE   = sums_e(end);

    % Subset cumulative sums at indices floor(n p)
    lorenz = sums_e(floor_np) ./ totalE;  % normalized L(p) in [0,1]
end


