# iszero on v::SparseVector calls zero(eltype(v))
function Base.iszero(a::AlgebraElement{M, <:AbstractMatrix}) where M
    v = coeffs(a)
    for idx in __nzind(v)
        iszero(v[idx]) || return false
    end
    return true
end

function Base.:(==)(
    X::AlgebraElement{A, <:AbstractMatrix},
    Y::AlgebraElement{A, <:AbstractMatrix},
) where A
    parent(X) === parent(Y) || return false

    suppx = __nzind(coeffs(X))
    suppy = __nzind(coeffs(Y))
    x_and_y = intersect(suppx, suppy)
    x_only = setdiff(suppx, suppy)
    for idx in union(suppx, suppy)
        if idx in x_and_y
            X[idx] == Y[idx] || return false
        elseif idx in x_only
            iszero(X[idx]) || return false
        else # idx in y_only
            iszero(Y[idx]) || return false
        end
    end
    return true
end

function add!(
    res::AlgebraElement{A, <:AbstractMatrix},
    X::AlgebraElement,
    Y::AlgebraElement
) where A
    suppx = __nzind(coeffs(X))
    suppy = __nzind(coeffs(Y))
    x_and_y = intersect(suppx, suppy)
    x_only = setdiff(suppx, suppy)
    for idx in union(suppx, suppy)
        if idx in x_and_y
            res[idx] = X[idx] + Y[idx]
        elseif idx in x_only
            res[idx] = X[idx]
        else # idx in y_only
            res[idx] = Y[idx]
        end
    end
    return res
end

function sub!(
    res::AlgebraElement{A, <:AbstractMatrix},
    X::AlgebraElement,
    Y::AlgebraElement
) where A
    suppx = __nzind(coeffs(X))
    suppy = __nzind(coeffs(Y))
    x_and_y = intersect(suppx, suppy)
    x_only = setdiff(suppx, suppy)
    for idx in union(suppx, suppy)
        if idx in x_and_y
            res[idx] = X[idx] - Y[idx]
        elseif idx in x_only
            res[idx] = X[idx]
        else # idx in y_only
            res[idx] = - Y[idx]
        end
    end
    return res
end

Base.:*(a::AbstractMatrix, X::AlgebraElement) = mul!(similar(X, a*X._elt), a, X)
Base.:*(X::AlgebraElement, a::AbstractMatrix) = mul!(similar(X, a*X._elt), X, a)
Base.:(/)(X::AlgebraElement, a::AbstractMatrix) = X * inv(a)

function mul!(
    res::AlgebraElement{A, <:AbstractMatrix},
    X::AlgebraElement,
    a::Union{Number, AbstractMatrix},
) where A
    for idx in __nzind(coeffs(X))
        res[idx] = X[idx]*a
    end
    return res
end

function mul!(
    res::AlgebraElement{A, <:AbstractMatrix},
    a::AbstractMatrix,
    X::AlgebraElement,
) where A
    for idx in __nzind(coeffs(X))
        res[idx] = a*X[idx]
    end
    return res
end

function fmac!(
    res::AbstractVector{<:AbstractMatrix},
    X::AbstractVector,
    Y::AbstractVector,
    mstr::MultiplicativeStructure,
)
    @assert res !== X
    @assert res !== Y
    for (j, y) in _nzpairs(Y)
        for (i, x) in _nzpairs(X)
            k = mstr[i, j]
            if k in __nzind(res)
                res[k] += x * y
            else
                res[k] = x*y
            end
        end
    end
    return res
end
