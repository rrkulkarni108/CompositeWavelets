%% Barbara: crossover plot for W1W2 vs single base
close all; clear; clc;
colormap gray; rng(0);

% User defined parameters 
M = 200;                         % number of simulations
sigma_vals = 5:2:40;             % noise std sweep 
show_images_for_one_sigma = false; % set true to visualize one case

%  load in image
A = double(imread('barbara_gray.bmp'));  % 512x512
[n, m] = size(A);  N = n*m;  peak = 255; 
fprintf('Image: %dx%d, N=%d\n', n, m, N);

% Filters 
hfilt_daub3 = [0.332670552950165   0.806891509311400   0.459877502118228 ...
               -0.135011020010067  -0.085441273882042   0.035226291882017];
hfilt_symm4 = [-0.075765714789273, -0.029635527645999,  0.497618667632015, ...
                0.803738751805916,  0.297857795605605, -0.099219543576935, ...
               -0.012603967262038,  0.032223100604043];

hfilt_haar = [sqrt(2)/2, sqrt(2)/2]; % Haar filter
hfilt_daub5 = [
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

hfilt_coif1 = [	.038580777748	-.126969125396	-.077161555496	...
				.607491641386	.745687558934	.226584265197	]; % coiflet 6 tap
hfilt_coif3  =  [	-.003793512864	.007782596426	.023452696142	...
				-.065771911281	-.061123390003	.405176902410	...
				.793777222626	.428483476378	-.071799821619	...
				-.082301927106	.034555027573	.015880544864	...
				-.009007976137	-.002574517688	.001117518771	...
				.000466216960	-.000070983303	-.000034599773	]; % coiflet 18 tap 
hfilt_coif4 = [	.000892313668	-.001629492013	-.007346166328	...
				.016068943964	.026682300156	-.081266699680	...
				-.056077313316	.415308407030	.782238930920	...
				.434386056491	-.066627474263	-.096220442034	...
				.039334427123	.025082261845	-.015211731527	...
				-.005658286686	.003751436157	.001266561929	...
				-.000589020757	-.000259974552	.000062339034	...
				.000031229876	-.000003259680	-.000001784985	]; % Coif 24 tap

%  Wavelet matrices 
W1   = Wavmat(hfilt_symm4, n, 2, 2);   % single base (Symmlet 4)
W2   = Wavmat(hfilt_coif3, n, 2, 2);   % second base (Daubechies 5)
W1W2 = W1 * W2;                        % product base

% Storage of values
S = numel(sigma_vals);
mse_single_mu = zeros(S,1);   mse_single_se = zeros(S,1);
mse_prod_mu   = zeros(S,1);   mse_prod_se   = zeros(S,1);
in_snr_db      = zeros(S,1);

%  loop over sigma 
for s = 1:S
    sigma = sigma_vals(s);
    tau = sigma * sqrt(2*log(N));   % universal threshold 
    in_snr_db(s) = 20*log10(std(A(:))/sigma);

    mse_single = zeros(M,1);
    mse_prod   = zeros(M,1);

    for mrep = 1:M
        Y = A + sigma * randn(n,m);

        % Single base (W1) 
        C1 = W1 * Y * W1';
        C1_th = C1 .* (abs(C1) > tau);
        A1 = W1' * C1_th * W1;
        mse_single(mrep) = mean((A1(:) - A(:)).^2);

        % Product base (W1W2) 
        C12 = W1W2 * Y * W1W2';
        C12_th = C12 .* (abs(C12) > tau);
        A12 = W1W2' * C12_th * W1W2;
        mse_prod(mrep) = mean((A12(:) - A(:)).^2);

        % visualize one case if desired
        if show_images_for_one_sigma && s==1 && mrep==1
            figure; 
            subplot(1,3,1); imagesc(Y); axis image off; title(sprintf('Noisy (\\sigma=%g)',sigma));
            subplot(1,3,2); imagesc(A1); axis image off; title('Single base W1');
            subplot(1,3,3); imagesc(A12); axis image off; title('Product base W1W2');
            drawnow;
        end
    end

    % Averages and standard errors
    mse_single_mu(s) = mean(mse_single);
    mse_single_se(s) = std(mse_single)/sqrt(M);
    mse_prod_mu(s)   = mean(mse_prod);
    mse_prod_se(s)   = std(mse_prod)/sqrt(M);

end

%  crossover detection = first sigma where product beats single
diff_mu = mse_prod_mu - mse_single_mu;  % negative means product better
idx = find(diff_mu < 0, 1, 'first');
if ~isempty(idx)
    sigma_star = sigma_vals(idx);
    snr_star = in_snr_db(idx);
    fprintf('Crossover at sigma* = %.3g (input SNR ≈ %.2f dB): product < single.\n', sigma_star, snr_star);
else
    sigma_star = NaN; snr_star = NaN;
    fprintf('No crossover in the tested sigma range.\n');
end

% plotting
%figure('Name','MSE vs sigma'); hold on; box on;
fig1 = figure('Name','MSE vs sigma'); hold on; box on;
errorbar(sigma_vals, mse_single_mu, 1.96*mse_single_se, '-o', 'DisplayName','Single (Symm4)');
errorbar(sigma_vals, mse_prod_mu,   1.96*mse_prod_se,   '-s', 'DisplayName','Product (Symm4 ∘ Coif3)');
yl = yline(min([mse_single_mu; mse_prod_mu]),':','Color',[.5 .5 .5]);
yl.Annotation.LegendInformation.IconDisplayStyle = 'off';



if ~isnan(sigma_star), 
    xl1 = xline(sigma_star,'--','\sigma^*','LabelVerticalAlignment','bottom'); 
    xl1.Annotation.LegendInformation.IconDisplayStyle = 'off'; end
xlabel('\sigma (noise std)'); ylabel('MSE'); title('Denoising MSE vs \sigma (Barbara, universal threshold)');
legend('Location','northwest');
exportgraphics(fig1, 'mse_vs_sigma.pdf', 'ContentType', 'vector');

%figure('Name','ΔMSE vs sigma'); hold on; box on;
fig2 = figure('Name','ΔMSE vs sigma'); hold on; box on;
plot(sigma_vals, diff_mu, '-d','DisplayName','MSE(product) - MSE(single)');
yl2 = yline(0,'k-'); 
yl2.Annotation.LegendInformation.IconDisplayStyle = 'off';


if ~isnan(sigma_star),
    xl = xline(sigma_star,'--','\sigma^*','LabelVerticalAlignment','bottom');
    xl.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
xlabel('\sigma'); ylabel('\Delta MSE'); title('Crossover plot (negative favors product)');
legend('Location','best');

exportgraphics(fig2, 'delta_mse_vs_sigma.pdf', 'ContentType', 'vector');



%  table of results
T = table(sigma_vals(:), in_snr_db(:), mse_single_mu, mse_prod_mu, diff_mu, ...
    'VariableNames', {'sigma','InputSNR_dB','MSE_Single','MSE_Product','DeltaMSE'});
disp(T);
