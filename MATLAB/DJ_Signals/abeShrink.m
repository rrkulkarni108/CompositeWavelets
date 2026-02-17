%% Helper function - ABE method (Figueiredo & Nowak)

function wd_abe = abeShrink(wd, sigma2)
% Amplitude-Scale Invariant Bayes Estimator 
% wd      : coefficient vector
% sigma2  : noise variance (use 1 if standardized)

wd = wd(:);
wd2 = wd.^2;

wd_abe = zeros(size(wd));
mask = wd2 > 3*sigma2;

wd_abe(mask) = (wd2(mask) - 3*sigma2) ./ wd(mask);
end
