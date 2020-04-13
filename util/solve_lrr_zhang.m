function [Z,E] = solve_lrr(X,A,X2,A2,lambda,lam,reg,alm_type,display)
% Aug 2013
% This routine solves the following nuclear-norm optimization problem,
% min |Z|_*+lambda*|E|_L
% s.t., X = AZ+E
% inputs:
%        X -- D*N data matrix, D is the data dimension, and N is the number
%             of data vectors.
%        A -- D*M matrix of a dictionary, M is the size of the dictionary
%        lambda -- parameter
%        reg -- the norm chosen for characterizing E, 
%            -- reg=0 (default),                  use the l21-norm 
%            -- reg=1 (or ther values except 0),  use the l1-norm
%        alm_type -- 0 (default)   use the exact ALM strategy
%                 -- 1             use the inexact ALM strategy
%               
if nargin < 9 || isempty(display)
    display = false;
end
if nargin<8 || isempty(alm_type)
    alm_type = 0 ;
end

if nargin<7 || isempty(reg)
    reg = 0;
end

Q = orth(A');
B = A*Q;

Q2 = orth(A2');
B2 = A2*Q2;

if reg==0
    if alm_type ==0 
        [Z,E] = exact_alm_lrr_l21v2(X,B,lambda,[],[],display);
    else
        [Z,E] = inexact_alm_lrr_l21(X,B,lambda,display);
    end
else
    if alm_type == 0
        [Z,E] = exact_alm_lrr_l1v2(X,B,lambda,[],[],display);
    else
        [Z,E] = inexact_alm_lrr_l2l2(X,B,X2,B2,lambda,lam,display);
    end
end
Z = Q*Z;