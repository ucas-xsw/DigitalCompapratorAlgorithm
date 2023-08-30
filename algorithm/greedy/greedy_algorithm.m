function [opt, sol] = greedy_algorithm(values, weights, capacity)
    n = length(values);
    ratio = values ./ weights;
    [~, indices] = sort(ratio, 'descend');
    sol = zeros(1, n);
    w = 0;
    v = 0;
    for i = 1:n
        if w + weights(indices(i)) <= capacity
            sol(indices(i)) = 1;
            w = w + weights(indices(i));
            v = v + values(indices(i));
        end
    end
    opt = v;
end
