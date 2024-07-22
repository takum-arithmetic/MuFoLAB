"""
    infnorm(x) = ||x||_∞ = maximum(abs.(x)) 

    @param x: vector
    @return maximum(abs.(x))
"""
function infnorm(x)
	m = typeof(x[1])(0.0)
	for i in x
		j = i < 0 ? -i : i
		m = j > m ? j : m
	end
	return m
end
