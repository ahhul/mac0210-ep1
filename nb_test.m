f = @(x) (x .^ 4) .- 1
df = @(x) 4 .* (x .^ 3)
newton_basins (f, df, 20, 20, 20);
