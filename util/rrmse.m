function r = rmse(H, Hest)
[P, K] = size(H);
error_all = zeros(P);
for i = 1:P
    for j = 1:P
%         error_all(i, j) = norm(H(i, :) - Hest(j, :))/sqrt(K);
        error_all(i, j) = sqrt((norm(H(i, :) - Hest(j, :),'fro').^2)/K);
    end
end
error_min = min(error_all, [], 2);
r = mean(error_min);