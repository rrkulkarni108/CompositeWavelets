function W = Wavmateps(h, N, eps, shift)
% Wavmat_eps -- Transform Matrix for ε-decimated Orthogonal DWT (periodic)
% Usage:
%   W = Wavmat_eps(h, N, eps, shift)
%
% Inputs
%   h      : orthonormal low-pass filter (row or column)
%   N      : signal length (power of 2)
%   shift  : integer alignment parameter like in Wavmat (Default: 2)
%   eps    : 1 x k0 vector with entries in {0,1}:
%            eps(j)=1 -> "odd" phase (same as standard DWT),
%            eps(j)=0 -> "even" phase (phase-flipped decimation)
%
% Output
%   W      : N x N orthogonal transform matrix (W' = W^{-1})
%
% Notes
%  * If you pass eps = ones(1,k0), you recover your original Wavmat layout.
%  * Inverse is W' (transpose), since construction is orthonormal.
%  * Follows the block-recursive construction from Vidakovic (1999), pp. 115–116.
%
% Example
%   dat = [1 0 -3 2 1 0 1 2]';
%   h   = [sqrt(2)/2 sqrt(2)/2];
%   W1  = Wavmat_eps(h, 8, [1 1 1], 2)  % standard DWT
%   W2  = Wavmat_eps(h, 8, [0 1 0], 2);   % ε-decimated variant
%   wt1 = W1*dat;  rec1 = W1'*wt1;  max(abs(rec1 - dat))
%   wt2 = W2*dat;  rec2 = W2'*wt2;  max(abs(rec2 - dat))
%
% Copyright (c) 2025 Brani+Assistant

    if nargin < 4 || isempty(shift)
        shift = 2;
    end
   k0=length(eps);
   J=log2(N);
    eps = mod(eps(:).', 2);   % force 0/1 row

    % --- Make QMF high-pass from low-pass (same convention as your Wavmat)
    h = h(:)'; 
    g = fliplr(conj(h) .* (-1).^(1:length(h)));

    % Zero-extend to length N for modular indexing (as in Wavmat)
    h = [h, zeros(1, N)];
    g = [g, zeros(1, N)];

    % Initialization
    oldmat = eye(2^(J - k0));  % identity in the coarsest space

    % Descend levels: k = k0, k0-1, ..., 1
    for k = k0:-1:1
        ubJk  = 2^(J - k);        % number of columns in the new h/g blocks
        ubJk1 = 2^(J - k + 1);    % number of rows in the new h/g blocks

        % --- Per-level phase control via an effective shift
        %     eps(k)=1 -> "odd" decimation (matches standard DWT)
        %     eps(k)=0 -> "even" decimation (phase-flipped)
        %
        %     We implement this by a per-level shift:
        %       shift_k = shift + (1 - eps(k));
        %     If your convention ends up flipped by 1, simply change to:
        %       shift_k = shift + eps(k);
        shift_k = shift + (1 - eps(k));

        % Build hmat, gmat (as in Wavmat), but with per-level shift_k
        hmat = zeros(ubJk1, ubJk);
        gmat = zeros(ubJk1, ubJk);
        for jj = 1:ubJk
            for ii = 1:ubJk1
                modulus = mod(N + ii - 2*jj + shift_k, ubJk1);
                if modulus == 0, modulus = ubJk1; end
                hmat(ii, jj) = h(modulus);
                gmat(ii, jj) = g(modulus);
            end
        end

        % Assemble next-level transform
        W = [oldmat * hmat'; gmat'];
        oldmat = W;
    end
    % Final W
    % size(W) == [N N]
end
