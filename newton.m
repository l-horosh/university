% Function for Newton Method
function res = newton(grad, hess, x0)
    epsilon = 0.001;
    alpha = 0.1;
    n = length(x0);
    prev_x = x0;
    x = x0 + 1;
    maxn = 1000;
    res = zeros(n, maxn);
    res(:, 1) = x0;
    k = 1;
    while k <= maxn && norm(prev_x - x) > epsilon
        prev_x = x; 
        x_opt = iteration(hess(prev_x),(-grad(x)),prev_x) + prev_x;
        x = prev_x + alpha*(x_opt - prev_x); % u_k+1 = u_k + alpha_k(u_k* - u_k)
        res(:, k) = x;
        k = k + 1;
    end
    res(:, k:end) = [];
    res = fliplr(res);
    res = res(:, 1);
    disp(['Number of iterations: ', num2str(n - 1)]);
end














% f3 = @(x) sin(x(1)) + sin(x(2));
% f_cos3 = @(x,y) sin(x) + sin(y);
% grad3 = @(x) [cos(x(1)); cos(x(2))];
% hess3 = @(x) [-sin(x(1)), 0; 0, -sin(x(2))];