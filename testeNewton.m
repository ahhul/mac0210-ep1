f = @(x) (x .^ 4) .- 1
df = @(x) 4 .* (x .^ 3)

newton (f, df, 100);
