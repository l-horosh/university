function [u_min, J_min] = cond_grad(rho, u0, eps, J, grad)
    % [val; point] = rho(l) - support function of U
    % A, c - parameters of quadratic minimization
    % J, grad - parameters in jeneral case 
    % u0 - initial approximation
    % eps - J tolerance
    n_max = 1000;
    u = u0;
    u_prev = u0;
    n = 0;
    % итерационная последовательность u_k+1 = uk + a_k(u_k* - u)
    while n == 0 || (abs(J(u) - J(u_prev)) >= eps && n < n_max)
        u_prev = u;
        [~, v] = rho(-grad(u_prev));
        delta = v - u; %u_k - u
        if prod(delta == 0)
            break;
        end
        alpha = fminbnd(@(x) J(u_prev + x * delta), 0, 1);
        u = u_prev + alpha * delta;
        n = n + 1;
    end
    u_min = u;
    J_min = J(u);
    if n == n_max
        disp("Algorithm diverged");
        u_min = [];
        J_min = [];
        return;
    else
        disp(['Algorithm converged in ', num2str(n), ' iteration(s)']);
        return;
    end   
end