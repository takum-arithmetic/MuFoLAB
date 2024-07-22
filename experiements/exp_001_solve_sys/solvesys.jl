"""
    Experiment 1: Compare solutions for different types of numbers for 
    solving linear systems of equations, Ax = b. In this experiment, 
    we compare the infinity norm of the difference between the true 
    solution and the computed solution for different types of numbers. 
    The types of numbers we compare are Takums, SoftPosits, BFloat16s, and 
    Floats.  We use the SuiteSparse Matrix Collection to test the different 
    number types.  The integer ids are the 'official' ident numbers 
    assigned by the collection. Currently id ∈ 1:3000.  The matrices are 
    cast to FULL.
"""

# Local dependencies
include("configs.jl")
include("xgesol.jl")

# Run the experiment
@timev f(x) = xgesol.(x,Ts)
@timev precomp = f.([905,906])
@timev results = f.(Ms) # map(f,Ms) better

# Plot and save results
matrix_data = hcat(results...);
df = DataFrame(matrix_data', string.(Ts))

avgs = [mean(df[:, i]) for i in 1:ncol(df)]';
df_stats = DataFrame(avgs, string.(Ts))

p = plot(title = "Sovling Ax = b", 
    xlabel = "Matrix #", 
    ylabel = "infinity Norm", 
    legend = :topleft, 
    size = (800, 600), 
    dpi = 300)

for i in 1:ncol(df)
    plot!(p, 1:nrow(df), df[:,i], label = names(df)[i], marker = :auto)
end
display(p)
savefig("results.png")
CSV.write("results.csv", df)