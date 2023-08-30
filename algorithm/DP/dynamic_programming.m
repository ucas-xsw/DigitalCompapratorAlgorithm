function [opt, sol] = dynamic_programming(values, weights, capacity)
    n = length(values);
    dp = zeros(n+1, capacity+1);
    for i = 1:n
        for j = 1:capacity+1
            if j >= weights(i)
                dp(i+1, j) = max(dp(i, j), dp(i, j-weights(i)) + values(i));
            else
                dp(i+1, j) = dp(i, j);
            end
        end
    end
    opt = dp(n+1, capacity+1);
    sol = zeros(1, n);
    j = capacity+1;
    for i = n:-1:1
        if dp(i+1, j) > dp(i, j)
            sol(i) = 1;
            j = j - weights(i);
        end
    end
end
