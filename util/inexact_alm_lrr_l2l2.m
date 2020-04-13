function [Z,E] = inexact_alm_lrr_l1(X,A,Y,B,lambda,lam,display)
% This routine uses Inexact ALM algorithm to solve the following nuclear-norm optimization problem:
% min |Z|_*+lambda*|E|_1
% s.t., X = AZ+E
% inputs:
%        X -- D*N data matrix, D is the data dimension, and N is the number
%             of data vectors.
%        A -- D*M matrix of a dictionary, M is the size of the dictionary

if nargin<4
    display = false;
end

tol = 1e-8;
maxIter = 1e6;
[d n] = size(X);
[d2 n2] = size(Y);
m = size(A,2);
m2 = size(B,2);
rho = 1.1;
max_mu = 1e10;
mu = 1e-6;
atx = A'*X;
aty = B'*Y;
inv_a = inv(A'*A+eye(m));
inv_b = inv(B'*B+eye(m2));
inv_ab = inv(A'*A+B'*B+eye(m));
%% Initializing optimization variables
% intialize
J = zeros(m,n);
Z = zeros(m,n);
E = sparse(d,n);
E2 = sparse(d2,n2);

Y1 = zeros(d,n);
Y2 = zeros(m,n);
Y3 = zeros(d1,n1);
%% Start main loop
iter = 0;
if display
    disp(['initial,rank=' num2str(rank(Z))]);
end
while iter<maxIter
    iter = iter + 1;
    %update J
    temp = Z + Y2/mu;
    [U,sigma,V] = svd(temp,'econ');
    sigma = diag(sigma);
    svp = length(find(sigma>1/mu));
    if svp>=1
        sigma = sigma(1:svp)-1/mu;
    else
        svp = 1;
        sigma = 0;
    end
    J = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
    %udpate Z
%     Z = inv_a*(atx-A'*E+J+(A'*Y1-Y2)/mu);
    Z = inv_ab*(atx-A'*E+aty-B'*E2+J+(A'*Y1-Y2)/mu+(B'*Y3-Y2)/mu);
    %update E
    xmaz = X-A*Z;
    temp = xmaz+Y1/mu;
%     E = max(0,temp - lambda/mu)+min(0,temp + lambda/mu);
    E = (temp/(2*lambda/mu+1));
    %update E2
    ymaz = Y-B*Z;
    temp2 = ymaz+Y3/mu;
%     E = max(0,temp - lambda/mu)+min(0,temp + lambda/mu);
    E2 = (temp2/(2*lam/mu+1));
    
    leq1 = xmaz-E;
    leq2 = Z-J;
    leq3 = ymaz-E2;
    stopC = max(max(max(abs(leq1))),max(max(abs(leq2))));
    stopC = max(stopC, max(max(abs(leq3))))
    if display && (iter==1 || mod(iter,50)==0 || stopC<tol)
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e') ...
            ',rank=' num2str(rank(Z,1e-3*norm(Z,2))) ',stopALM=' num2str(stopC,'%2.3e')]);
    end
    if stopC<tol 
        break;
    else
        Y1 = Y1 + mu*leq1;
        Y2 = Y2 + mu*leq2;
        Y3 = Y3 + mu*leq3;
        mu = min(max_mu,mu*rho);
    end
end
