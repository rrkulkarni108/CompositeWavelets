# Composite Wavelet Matrix Based Transforms and Applications 
Radhika Kulkarni and Brani Vidakovic

This repository contains MATLAB and Python files of simulations and examples in the paper. Figures and simulations can be reproduced using these files. 
Please reach out to rrkulkarni@tamu.edu with any inquiries. 


### Abstract

Orthogonal wavelet transforms are a cornerstone of modern signal and image denoising
because they combine multiscale representation, energy preservation, and perfect reconstruction. In this
paper, we show that these advantages can be retained and substantially enhanced by moving beyond classical
single–basis wavelet filterbanks to a broader class of composite wavelet–like matrices. By combining
orthogonal wavelet matrices through products, Kronecker products, and block–diagonal constructions, we
obtain new unitary transforms that generally fall outside the strict wavelet filterbank class, yet remain fully
invertible and numerically stable.

The central finding is that such composite transforms induce stronger concentration of signal energy
into fewer coefficients than conventional wavelets. This increased sparsity, quantified using Lorenz curve
diagnostics, directly translates into improved denoising under identical thresholding rules. Extensive sim-
ulations on Donoho–Johnstone benchmark signals, complex–valued unitary examples, and adaptive block
constructions demonstrate consistent reductions in mean–squared error relative to single–basis transforms.
Applications to atmospheric turbulence measurements and image denoising of the Barbara benchmark
further confirm that composite transforms better preserve salient structures while suppressing noise.

From an engineering perspective, the proposed framework is attractive because it requires no modification
to the standard transform–shrink–inverse pipeline and preserves orthogonality and perfect reconstruction
by construction. The results suggest a shift in viewpoint: wavelet transforms can be treated as modular
orthogonal operators that may be algebraically combined to better match signal structure. Composite
wavelet–like matrices therefore provide a practical and flexible extension of classical wavelet methods for
denoising, compression, and multiscale analysis in one– and higher–dimensional settings
