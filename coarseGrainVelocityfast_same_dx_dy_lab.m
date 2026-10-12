function [Vx_grid, Vy_grid] = coarseGrainVelocity( ...
    x, y, vx, vy, grid_size_x, grid_size_y, sigma, X, Y)

Vx_grid = zeros(grid_size_y, grid_size_x);
Vy_grid = zeros(grid_size_y, grid_size_x);

Np = length(x);

for i = 1:grid_size_y
    for j = 1:grid_size_x

        w_sum = 0;
        vx_sum = 0;
        vy_sum = 0;

        for p = 1:Np

            dx = X(i,j) - x(p);
            dy = Y(i,j) - y(p);

            r2 = dx^2 + dy^2;

            % LOCALITY CONDITION
            if r2 < (3*sigma)^2

                w = exp(-r2/(2*sigma^2));

                vx_sum = vx_sum + vx(p)*w;
                vy_sum = vy_sum + vy(p)*w;
                w_sum = w_sum + w;

            end
        end

        if w_sum > 1e-8
            Vx_grid(i,j) = vx_sum / w_sum;
            Vy_grid(i,j) = vy_sum / w_sum;
        end

    end
end

end