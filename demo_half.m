% load('data.mat')
addpath('./util/')

% define random states
% rand('state',1);
% randn('state',1);

L = size(Y,1);
N = size(Y,2);
nEnd = size(M,2);

p.sources = nEnd;
% p.derta = 0;
%p.lambda = 0.06;
% p.lambda = 1./sqrt(L)*sum((sqrt(N)-sum(Y,2)./sum(Y.^2,2))./sqrt(N-1));
% 
% [Wvca, location, y] = VCA(Y, 'Endmembers', nEnd);
% Winit = Wvca;
% Hinit = max(pinv(Winit) * Y, 1e-10);
% Hinit = Hinit./repmat(sum(Hinit), nEnd, 1);
p.S = Hinit;

% p.lambda=0;
% dertavalue = [5,15];
% for ii=1:length(dertavalue)
%     p.derta = dertavalue(ii);
%     [Xest, Mest] = nnschalfprior3(Y, Winit, p, 'result');
%     Sam = sam(M, Mest); 
%     r= rmse(X, Xest, Sam(1,:), Sam(2,:));
%     aa(ii) = Sam(3,nEnd+1);
%     rr(ii) = r(3,nEnd+1);
% %     figure;
% %     for i=1:nEnd
% %         plot(Mest(:,i));hold on;
% %     end
% end

p.derta = 5;
% lambdavalue = [5e-4,1e-3,5e-3,0.01,0.05,0.1,0.2,0.5,1];
% lambdavalue = [0.01,0.1,0.5,1];
% lambdavalue = [10,13,15];
% lambdavalue = [0.3,0.4,0.]6;
lambdavalue = [0.1];
for ii=1:length(lambdavalue)
    p.lambda = lambdavalue(ii);
    [Xest, Mest] = nnschalfprior3(Y, Winit, p, 'result');
    fprintf('initial SAD estimation:\n');
    Sam = sam(M, Mest); 
    fprintf('initial RMSE estimation:\n');
    r = rmse(X, Xest, Sam(1,:), Sam(2,:));
    aa(ii) = Sam(3,nEnd+1);
    rr(ii) = r(3,nEnd+1);
end


% [Xest, Mest] = nnschalfprior3(Y, Winit, p, 'result');
% angle = sadd(M, Mest)
% r = rmse(X, Xest)
%    
figure;
for i=1:nEnd
    subplot(2,nEnd,i);imagesc(reshape(Xest(i,:),nc,nl));
    subplot(2,nEnd,i+nEnd);imagesc(reshape(X(i,:),nc,nl));
end

figure;
for i=1:nEnd
    plot(Mest(:,i));hold on;
end

figure;
for i=1:nEnd
    plot(M(:,i));hold on;
end