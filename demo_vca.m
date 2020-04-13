% load('data.mat')
addpath('./util/')

% define random states
rand('state',10);
randn('state',10);

L = size(Y,1);
N = size(Y,2);
nEnd = size(M,2);

% for ii=1:20
%     [Wvca, location, y] = VCA(Y, 'Endmembers', nEnd);
%     Winit = Wvca;
%     Hinit = max(pinv(Winit) * Y, 1e-10);
%     Hinit = Hinit./repmat(sum(Hinit), nEnd, 1);
% 
%     Sam = sam(M, Winit); 
%     r = rmse(X, Hinit, Sam(1,:), Sam(2,:));
%     aa(ii) = Sam(3,nEnd+1);
%     rr(ii) = r(3,nEnd+1);
% end
% mean(aa)
% mean(rr)

[Wvca, location, y] = VCA(Y, 'Endmembers', nEnd);
Winit = Wvca;
Hinit = max(pinv(Winit) * Y, 1e-10);
Hinit = Hinit./repmat(sum(Hinit), nEnd, 1);

Sam = sam(M, Winit); 
rmse(X, Hinit, Sam(1,:), Sam(2,:));
% 
% % 
% figure(1);
% for i=1:nEnd
%     subplot(2,nEnd,i);imagesc(reshape(Hinit(i,:),nc,nl));
%     subplot(2,nEnd,i+nEnd);imagesc(reshape(X(i,:),nc,nl));
% end
% 
% figure(2)
% for i=1:nEnd
%     plot(Winit(:,i));hold on;
% end
% 
% figure(3)
% for i=1:nEnd
%     plot(M(:,i));hold on;
% end