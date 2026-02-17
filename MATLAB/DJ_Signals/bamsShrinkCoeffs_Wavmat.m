%% Helper function - BAMS levelwise
% Based on demo2d, BAMS.m, Wavmat.m and Angelini Vidakovic gamma 

function [wd_shrunk, pars] = bamsShrinkCoeffs_Wavmat(wd, n, nl)
% BAMS/Angelini shrinkage applied to coefficient vector from Wavmat ordering.
% Uses bayesrule.m (Vidakovic & Ruggeri).
%
% wd: n x 1 coeff vector (for ex. wd = W*s')
% n : length (power of 2)
% nl: depth used in Wavmat(h,n,nl,shift)

J = log2(n);
if abs(J-round(J)) > 1e-12, error('n must be power of 2'); end

wd_shrunk = wd;

% Wavmat depth nl => scaling block length = 2^(J-nl)
% Coarsest detail index in indexing = i = J-nl
coarsest = max(J - nl, 1);
finest   = J - 1;

%  estimate mu from finest detail block (last half)
indx = (2^(J-1)+1) : 2^J;
niz  = wd(indx);
nizo = sort(niz);
mm   = length(nizo);
low  = max(floor(0.25*mm),1);
high = max(floor(0.75*mm),1);
pseudos = abs(nizo(high) - nizo(low)) / 1.5;   % Tukey pseudo-sigma
mu = 1 / max(pseudos^2, 1e-8);

%  set tau and eps levelwise (same forms as BAMS.m)
tauu = sqrt(max(std(wd)^2 - 1/mu, 1e-4));

eps_vec = nan(1, finest);
tau_vec = nan(1, finest);

for i = finest:-1:coarsest
    eps_i = 1 - 1/(i - coarsest + 1)^1.5;
    tau_i = tauu*(i+1)^2;

    ind = (2^i + 1) : 2^(i+1);
    a = wd(ind);

    b = zeros(size(a));
    for j = 1:length(a)
        [~, b(j)] = bayesrule(a(j), mu, tau_i, eps_i);  % b(j) = bre
    end

    wd_shrunk(ind) = b;
    eps_vec(i) = eps_i;
    tau_vec(i) = tau_i;
end

pars.mu = mu; pars.tau = tau_vec; pars.eps = eps_vec;
pars.coarsest = coarsest; pars.finest = finest;

end
