function [res,angle_all] = angle(W, West)
P = size(W, 2);
Q = size(West, 2);
cor_all = zeros(P, Q);
for i = 1:P
    for j = 1:Q
        cor_all(i, j) = W(:, i)'*West(:, j)/(norm(W(:, i))*norm(West(:, j)));
    end
end
cor_max = max(cor_all, [], 2);
angle_all = acos(cor_max);
res = mean(angle_all);
end