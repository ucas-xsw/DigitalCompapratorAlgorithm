function [X] = nelders_favorite_initial_simplex(initial_point,v_diag)
% -- ret = nelders_favorite_func(x)
%
%     Computes initial simplex for Nelder's favorite function, as described in
%     [1], Excercise 8.6.1.
%
%   References:
%     [1] Kelley, Carl T. Iterative methods for optimization. Society for
%         Industrial and Applied Mathematics, 1999.

  X = [initial_point, diag(v_diag)];
end
