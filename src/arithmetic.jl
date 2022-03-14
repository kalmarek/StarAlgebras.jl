# module structure:
Base.:*(a::Number, X::AlgebraElement) =
    mul!(similar(X, a*X._elt), X, a)
Base.:*(X::AlgebraElement, a::Number) = a * X
Base.:(/)(X::AlgebraElement, a::Number) = inv(a) * X

# TODO: handle this through mul!?
Base.:(//)(X::AlgebraElement, a::Number) =
    AlgebraElement(coeffs(X) .// a, parent(X), X._elt//a)

# ring structure:
Base.:-(X::AlgebraElement) = neg!(similar(X), X)

function _preallocate_output(op, X::AlgebraElement, Y::AlgebraElement)
    _elt = op(X._elt, Y._elt)
    if _elt isa AbstractMatrix
        w = similar(X, _elt)
        _zero!(coeffs(w))
        return w
    else
        return coeffs(Y) isa DenseArray ? similar(Y, _elt) : similar(X, _elt)
    end
end

Base.:+(X::AlgebraElement, Y::AlgebraElement) = add!(_preallocate_output(+, X, Y), X, Y)
Base.:-(X::AlgebraElement, Y::AlgebraElement) = sub!(_preallocate_output(-, X, Y), X, Y)
Base.:*(X::AlgebraElement, Y::AlgebraElement) = mul!(_preallocate_output(*, X, Y), X, Y)

Base.:^(a::AlgebraElement, p::Integer) = Base.power_by_squaring(a, p)

# mutable API; TODO: replace with MutableArithmetic

__nzind(v::AbstractVector) = eachindex(v)
__nzind(v::AbstractSparseVector) = SparseArrays.nonzeroinds(v)

function _zero!(v::AbstractVector, elt)
    for idx in __nzind(v)
        v[idx] = zero(elt)
    end
    return v
end

zero!(a::AlgebraElement) = (_zero!(coeffs(a), a._elt); a)

function neg!(res::AlgebraElement, X::AlgebraElement)
    @assert parent(res) === parent(X)
    res.coeffs .= -coeffs(X)
    return res
end

function add!(res::AlgebraElement, X::AlgebraElement, Y::AlgebraElement)
    @assert parent(res) === parent(X) === parent(Y)
    # res = (res === X || res === Y) ? similar(res) : res
    res.coeffs .= coeffs(X) .+ coeffs(Y)
    return res
end

function sub!(res::AlgebraElement, X::AlgebraElement, Y::AlgebraElement)
    @assert parent(res) === parent(X) === parent(Y)
    # res = (res === X || res === Y) ? similar(res) : res
    res.coeffs .= coeffs(X) .- coeffs(Y)
    return res
end

function mul!(res::AlgebraElement, X::AlgebraElement, a::Number)
    @assert parent(res) === parent(X)
    res.coeffs .= a .* coeffs(X)
    return res
end

function mul!(
    res::AbstractVector,
    X::AbstractVector,
    Y::AbstractVector,
    ms::MultiplicativeStructure,
)
    res = (res === X || res === Y) ? zero(res) : _zero!(res, zero(eltype(res)))
    return fmac!(res, X, Y, ms)
end

function mul!(res::AlgebraElement, X::AlgebraElement, Y::AlgebraElement)
    @assert parent(res) === parent(X) === parent(Y)
    res = (res === X || res === Y) ? zero(res) : zero!(res)
    return fmac!(res, X, Y)
end

function fmac!(res::AlgebraElement, X::AlgebraElement, Y::AlgebraElement)
    fmac!(coeffs(res), coeffs(X), coeffs(Y), parent(res).mstructure)
    return res
end

_nzpairs(v::AbstractVector) = pairs(v)
_nzpairs(v::AbstractSparseVector) =
    zip(SparseArrays.nonzeroinds(v), SparseArrays.nonzeros(v))

function fmac!(
    res::AbstractVector,
    X::AbstractVector,
    Y::AbstractVector,
    mstr::MultiplicativeStructure,
)
    @assert res !== X
    @assert res !== Y
    for (j, y) in _nzpairs(Y)
        for (i, x) in _nzpairs(X)
            res[mstr[i, j]] += x * y
        end
    end
    return res
end
