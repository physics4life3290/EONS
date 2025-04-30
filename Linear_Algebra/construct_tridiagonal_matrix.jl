




# General Algorithm/Function
function build_tridiagonal(a::Vector, b::Vector, c::Vector)
    n = length(b)
    T = zeros(Float64, n, n)  # Create a zero matrix of size n x n
    
    for i in 1:n
        T[i, i] = b[i]  # Assign the diagonal element
        if i > 1
            T[i, i-1] = a[i]  # Assign the subdiagonal element
        end
        if i < n
            T[i, i+1] = c[i]  # Assign the superdiagonal element
        end
    end
    
    return T
end
