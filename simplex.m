% Function for Simplex Method
function [umin, Jmin] = simplex(A, b, c, u_0)
    Jmin = -inf;
    umin=[];
    r=rank(A);
    eps = 0.001;
    i = 0;
    while true
        i = i+1;
        [Jbase, ~] = find(u_0 >= eps);
        [Jnbase, ~] = find(u_0 < eps);
        if numel(Jbase) < r
            m = r - numel(Jbase);
            Jbase = [Jbase; Jnbase(1:m)];
            Jnbase = Jnbase(m+1:end);
        end
        B = A(1:end, Jbase);
        F = A(1:end, Jnbase);
        cb = c(Jbase);
        cf = c(Jnbase);
        x = u_0(Jbase);
        C = B \ F;
        delta = C' * cb - cf;
        Jnbpos = find(delta > eps); 
        if isempty(Jnbpos)
            break;
        end
        k = min(Jnbpos);
        gamma = C(1:end,k);
        help = find(gamma > eps);
        if isempty(help)
            umin = [];
            Jmin = -inf;
            return;
        end
        theta = min(x(help) ./ gamma(help));
        umin = zeros(size(u_0));
        umin(Jbase) = x - theta * gamma;
        umin(Jnbase(k)) = theta;
        Jmin = umin' * c;
        u_0 = umin;
    end
    disp(['Algorithm converged in ', num2str(i), ' iteration(s)']);
end