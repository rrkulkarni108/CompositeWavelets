%% Helper function - BAMS global method 

function [wd_shrunk, pars] = bamsShrinkGlobal(wd)
% Global BAMS shrinkage using bayesrule.m
% Uses one (mu,tau,eps) for all coefficients without any  blocks fpr the
% levels
wd = wd(:);
n = length(wd);

% setting the parameters: mu 
w = sort(wd);
low  = max(floor(0.25*n),1);
high = max(floor(0.75*n),1);
pseudos = abs(w(high) - w(low)) / 1.5; %tukey's book
mu = 1 / max(pseudos^2, 1e-8);

% setting pars: tau and epsilon (NOT levelwise)
tau = sqrt(max(std(wd)^2 - 1/mu, 0.0001));

% eps is fixed default like demo2d figure/ Angelini also uses ~0.95
eps = 0.95;

wd_shrunk = zeros(size(wd));
for j = 1:n
    [~, wd_shrunk(j)] = bayesrule(wd(j), mu, tau, eps);
end

pars.mu = mu; pars.tau = tau; pars.eps = eps;
end
