function [W, H] = nmfupk(Y, Winit, Hinit, prior, parameter)
% Nonnegative Matrix Factorization for Hyperspectral Unmixing Using Prior
% Knowledge of Spectral Signatures 
%
% [NMFupk] Wei Tang, Zhenwei Shi, and Zhenyu An, "Nonnegative matrix 
% factorization for hyperspectral unmixing using prior knowledge of 
% spectral signatures", Optical Engineering, 51(8), 087001 (August 2012).
%
% -------------------------------------------------------------------
% Usage:
%
% [W, H] = nmfupk(Y, Winit, Hinit, prior, parameter)
%
% ------- Input variables -------------------------------------------
%
%  Y - hyperspectral data matrix with dimensions L(bands) x K(pixels)
%
%  Winit - initial endmember matrix with dimensions L(bands) x 
%  P(endmember number)
%
%  Hinit - initial abundance matrix with dimensions P(endmember number) x 
%  K(pixels)
%
%  prior - the vector of locations of prior endmembers in Winit
%
%  parameter.
%            maxiter - maximum iteration number, 1 x 1 matrix, default: 200
%            tol - error tolerance, default: 1e-3
%            isSumToOne - whether the sum-to-one constraint is forced,
%            default: false
%
% ------- Output variables -------------------------------------------
%
% W - estimated endmember matrix: the prior endmembers are the left-most 
% columns
%
% H - estimated abundance matrix
%
% ---------------------------------------------------------------------
%
% Please see [NMFupk] for more details.
%
% Please contact Wei Tang (tangwei@sa.buaa.edu.cn) to report bugs or 
% provide suggestions and discussions for the codes.
%
% ---------------------------------------------------------------------
% version: 1.0 (30-Mar-2014)
% ---------------------------------------------------------------------
%
% Copyright (Mar, 2014):       Wei Tang (tangwei@sa.buaa.edu.cn)
%                              Zhenwei Shi (shizhenwei@buaa.edu.cn)
%                              Zhenyu An (an.zhenyu.buaa@gmail.com)
%
% NMFupk is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ---------------------------------------------------------------------

%% Check input
if isfield(parameter, 'isSumToOne')
    isSumToOne = parameter.isSumToOne;
else
    isSumToOne = false;
    warning('The sum-to-one constraint is excluded!');
end

if isfield(parameter, 'maxiter')
    maxiter = parameter.maxiter;
else
    maxiter = 200;
    warning('Using default maximum iteration: 200');
end

if isfield(parameter, 'tol')
    tol = parameter.tol;
else
    tol = 1e-3;
    warning('Using default error tolerance: 0.001');
end
[L, K] = size(Y);
P = size(Winit, 2);
if L ~= size(Winit, 1) || P ~= size(Hinit, 1) || K ~= size(Hinit, 2)
    error('The sizes of the input data are incompatible!');
end

%% Prepare data
noprior = setdiff([1:P], prior);
W1 = Winit(:, prior);
W2 = Winit(:, noprior);
H1 = Hinit(prior, :);
H2 = Hinit(noprior, :);

error_old = 0.5*norm(Y - W1*H1 - W2*H2, 'fro')^2;
epsilon = 1e-16;

W1tW1 = W1'*W1;
%% Main iteration
for iter = 1:maxiter
    
    W1(W1 < 0) = epsilon;
    W2(W2 < 0) = epsilon;
    H1(H1 < 0) = epsilon;
    H2(H2 < 0) = epsilon;
         
    if (rem(iter, 10) == 0) && isSumToOne
        sumH = sum([H1; H2]);
        H1 = H1./repmat(sumH, size(H1, 1), 1);
        H2 = H2./repmat(sumH, size(H2, 1), 1);
    end
    
    W2 = W2.*(Y*H2')./(W1*H1*H2'+W2*(H2*H2') + epsilon);
    H1 = H1.*(W1'*Y)./(W1tW1*H1+W1'*W2*H2 + epsilon);
    H2 = H2.*(W2'*Y)./(W2'*W1*H1+W2'*W2*H2 + epsilon);
    
    error_new = 0.5*norm(Y - W1*H1 - W2*H2, 'fro')^2;
    
    toltemp = abs(error_old - error_new)/error_old;
    
    error_old = error_new;
    
    if rem(iter, 10) == 0 || iter == 1
        display(['iter: ', num2str(iter), ' error: ', ...
            num2str(error_new), ' tol: ', num2str(toltemp)]);
    end
    
    if toltemp < tol
        break;
    end
end
if iter == maxiter
    display('Maximum iteration has been reached!');
end

W = [W1, W2];
H = [H1; H2];
end