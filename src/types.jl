abstract type AbstractStarAlgebra{O,T} end

struct StarAlgebra{O,T,M<:MultiplicativeStructure,B<:AbstractBasis{T}} <:
       AbstractStarAlgebra{O,T}
    object::O
    mstructure::M
    basis::B

    function StarAlgebra(
        obj,
        basis::AbstractBasis,
        mstr::MultiplicativeStructure,
    )
        O = typeof(obj)
        T = eltype(basis)
        M = typeof(mstr)
        B = typeof(basis)

        return new{O,T,M,B}(obj, mstr, basis)
    end

    function StarAlgebra(obj, mstr::MultiplicativeStructure)
        O = typeof(obj)
        T = Symbol
        M = typeof(mstr)
        B = Basis{T,Int}

        return new{O,T,M,B}(obj, mstr)
    end
end

# TrivialMStructure:
StarAlgebra(obj, basis::AbstractBasis) = StarAlgebra{false}(obj, basis)

function StarAlgebra{Tw}(obj, basis::AbstractBasis) where {Tw}
    mstr = TrivialMStructure{Tw}(basis)
    return StarAlgebra(obj, basis, mstr)
end

# CachedMStructure:
function StarAlgebra(
    obj,
    basis::AbstractBasis,
    cache_size::Tuple{<:Integer,Integer};
    precompute = false,
)
    return StarAlgebra{false}(obj, basis, cache_size, precompute = precompute)
end

function StarAlgebra{Tw}(
    obj,
    basis::AbstractBasis,
    cache_size::Tuple{<:Integer,Integer};
    precompute = false,
) where {Tw}
    mstr = CachedMTable{Tw}(basis, table_size = cache_size)
    precompute && complete!(mstr)
    return StarAlgebra(obj, basis, mstr)
end

hasbasis(A::StarAlgebra) = isdefined(A, :basis)

basis(A::StarAlgebra) = A.basis
object(A::StarAlgebra) = A.object
# Base.eltype(A::StarAlgebra{O,B}) where {O,B} = eltype(B)

struct AlgebraElement{A,T,V<:AbstractVector{T}}
    coeffs::V
    parent::A
    _elt::T

    function AlgebraElement(coeffs::AbstractVector, A::AbstractStarAlgebra, cf = first(coeffs))
        if hasbasis(A)
            @assert length(coeffs) == length(basis(A))
        end
        return new{typeof(A),eltype(coeffs),typeof(coeffs)}(coeffs, A, cf)
    end
end

coeffs(a::AlgebraElement) = a.coeffs
Base.parent(a::AlgebraElement) = a.parent
Base.eltype(a::AlgebraElement) = eltype(coeffs(a))

### constructing elements

function Base.zero(A::AbstractStarAlgebra, _elt = 1)
    if hasbasis(A)
        I = SparseArrays.indtype(basis(A))
        return AlgebraElement(sparsevec(I[], eltype(_elt)[], length(basis(A))), A, _elt)
    end
    return throw(
        "Algebra without basis; to construct zero use the `AlgebraElement` constructor directly.",
    )
end

function Base.one(A::AbstractStarAlgebra, _elt = 1)
    hasbasis(A) && return A(one(object(A)), _elt)
    return throw(
        "Algebra without basis; to construct one use the `AlgebraElement` constructor directly.",
    )
end

Base.zero(a::AlgebraElement) = (b = similar(a); return zero!(b))
Base.one(a::AlgebraElement) = one(parent(a), a._elt)
Base.iszero(a::AlgebraElement) = iszero(coeffs(a))

function Base.isone(a::AlgebraElement)
    b = basis(parent(a))
    k = findfirst(!iszero, coeffs(a))
    k === nothing && return false
    isone(a[k]) || return false
    return isone(b[k]) && findnext(!iszero, coeffs(a), k+1) === nothing
end

function (A::AbstractStarAlgebra{O,T})(elt::T, _elt = 1) where {O,T}
    if hasbasis(A)
        b = basis(A)
        i = b[elt]
        return AlgebraElement(sparsevec([i], [one(_elt)], length(b)), A, _elt)
    else
        throw("Algebra without basis: cannot coerce $elt.")
    end
end

function (A::AbstractStarAlgebra)(x::Number)
    g = one(object(A))
    res = A(g, x)
    res[g] *= x
    return res
end

Base.similar(X::AlgebraElement, _elt = X._elt) =
    AlgebraElement(similar(coeffs(X), typeof(_elt)), parent(X), _elt)

function AlgebraElement{T}(X::AlgebraElement) where {T}
    v = coeffs(X)
    w = similar(v, T)
    _elt = convert(T, X._elt)
    w .= v
    return AlgebraElement(w, parent(X), _elt)
end
