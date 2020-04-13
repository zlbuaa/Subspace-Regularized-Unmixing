% load('data.mat')
addpath('./util/')

% define random states
rand('state',1);
randn('state',1);

L = size(Y,1);
N = size(Y,2);
nEnd = size(M,2);

p.sources = nEnd;
p.derta = 5;
% p.lambda = 0.06;
% % p.lambda = 1./sqrt(L)*sum((sqrt(N)-sum(Y,2)./sum(Y.^2,2))./sqrt(N-1));
% p.mu = 1e-4;

[Wvca, location, y] = VCA(Y, 'Endmembers', nEnd);
Winit = Wvca;
Hinit = max(pinv(Winit) * Y, 1e-10);
Hinit = Hinit./repmat(sum(Hinit), nEnd, 1);
fprintf('initial SAD estimation:\n');
Sam = sam(M, Winit); 
fprintf('initial RMSE estimation:\n');
r = rmse(X, Hinit, Sam(1,:), Sam(2,:));
% [Z,E] = solve_lrr(Y,Y,0.01,1,1,1);

% lambdavalue = [5e-4,1e-3,5e-3,0,0.01,0.05,0.1,0.2,0.3];
% lambdavalue = [5e-4,1e-3,5e-3,0.01,0.05,0.1,0.2,0.3];
% lambdavalue = [0.01,0.1,0.3];
lambdavalue = 15;
% muvalue = [1e-4,1e-3,0.01,0.1,0.2,0.5];
muvalue = [0.5,1,1.5,5,10];
% muvalue = [5e-4,5e-3,5e-2];
% muvalue = 0.6;


for ii = 1:length(muvalue)
    for jj = 1:length(lambdavalue)
        p.lambda = lambdavalue(jj);
        p.mu = muvalue(ii);
        
        p.S = Hinit;

        ZZ=eye(nc*nl)-Z;
        [Sest, Aest] = newregu(Y, Winit, ZZ, p, 'result1');
        fprintf('initial SAD estimation:\n');
        Sam = sam(M, Aest); 
        fprintf('initial RMSE estimation:\n');
        r = rmse(X, Sest, Sam(1,:), Sam(2,:));

        aa(ii,jj) = Sam(3,nEnd+1);
        rr(ii,jj) = r(3,nEnd+1);
    end
end

