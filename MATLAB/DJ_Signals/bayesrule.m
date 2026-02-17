%% Helper function - bayesrule 

function [me, bre]=bayesrule(d, mu, tau, eps)
%
%
%  Bayesrule: Calculates Marginal Distribution and Shrinkage Bayes Rule
%  (needed by BAMS)
%  Usage
%   [me, br]= bayesrule(d, mu, tau, eps)
%  Inputs
%    d      Wavelet Coefficient
%    mu     Parameter - concerning prior on the scale of noise
%    tau    Scale of the prior on location of d (signal)
%    eps    Weight of point mass at zero in the prior on signal
%  Outputs
%    me     value of the marginal distribution at d
%    br     value of the bayes rule at d
%
de = 1/2 * sqrt(2 * mu) * exp(- abs(d)*sqrt(2 * mu) );

m = (tau * exp(-abs(d)/tau) - 1/sqrt(2 * mu) * ...
   exp( - sqrt(2 * mu) * abs(d) ) )/(2 * tau^2 - 1/mu);
me = eps * de + (1-eps) * m;

if d >= 0,
    br= ( tau * (tau^2 - 1/(2*mu)) * d * exp(-d/tau) + ...
    tau^2/mu * ( exp(- sqrt(2*mu)*d) - exp(-d/tau)))/  ...
    ( (tau^2 - 1/(2*mu)) * (tau *  exp(-d/tau) -  ...
    1/sqrt(2*mu)*  exp(- sqrt(2*mu)*d) ) ); 
else
    br= ( tau * (tau^2 - 1/(2*mu)) * d * exp(d/tau) - ...
    tau^2/mu * ( exp(sqrt(2*mu)*d) - exp(d/tau)))/  ...
    ( (tau^2 - 1/(2*mu)) * (tau *  exp(d/tau) -  ...
    1/sqrt(2*mu)*  exp(sqrt(2*mu)*d) ) ); 
end

bre= (1-eps)*m*br / ( (1-eps)*m + eps * de );

%
% Copyright (c) 2000  B. Vidakovic and F. Ruggeri
%


