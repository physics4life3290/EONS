




function build_spherical_poisson_matrix(n, Δr)
    # Allocate a, b, and c only once
    a = Vector{Float64}(undef, n)
    b = Vector{Float64}(undef, n)
    c = Vector{Float64}(undef, n)
    
    # Loop over all indices
    for i in 1:n
        if i == 1
            b[i] = -2/Δr^2
            c[i] = +2/Δr^2
        else
            r = i * Δr
            a[i] = 1/Δr^2 - 1/(r * Δr)   # sub-diagonal
            b[i] = -2/Δr^2              # diagonal
            c[i] = 1/Δr^2 + 1/(r * Δr)  # super-diagonal
        end
    end

    # Return the tridiagonal matrix along with the coefficients for the super-diagonal
    return build_tridiagonal(a, b, c), c
end