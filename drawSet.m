function drawSet(rho, N)
    i = 1 : N;
    L = [cos(2 * pi * i / N); sin(2 * pi * i / N)];
    vals = zeros(N, 1);
    points = zeros(2, N);
    for i = 1 : N
        [vals(i), points(:, i)] = rho(L(:, i));
    end
    plot([points(1, :), points(1, 1)], [points(2, :), points(2, 1)], 'black', 'LineWidth', 2.5);
end