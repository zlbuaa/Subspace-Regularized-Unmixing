function [loc_prior, matched_num] = match(Wp, Winit)
P = size(Winit, 2);%number of endmembers
prior_num = size(Wp, 2);%number of known endmembers

matched_num = 0;
loc_prior = zeros(1, prior_num);
cor_all = zeros(prior_num, P);
for i = 1:prior_num
    for j = 1:P    
        if isempty(find(loc_prior == j)) == 1
            cor_all(i, j) = Wp(:, i)'*Winit(:, j)...
                /(norm(Wp(:, i))*norm(Winit(:, j)));
        else
            cor_all(i, j) = 0;
        end
    end
    [val, loc] = max(cor_all(i, :));
    loc_prior(i) = loc;
    if val > 0.9
        matched_num = matched_num + 1;
    end
end

loc_prior = sort(loc_prior);
end