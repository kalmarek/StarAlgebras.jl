@testset "Algebra and Elements" begin
    A = [:a, :b, :c]
    b = StarAlgebras.Basis{UInt8}(words(A, radius=2))
    l = length(b)

    RG = StarAlgebras.StarAlgebra(one(first(b)), b, (4, 4))

    a = rand(l)

    @test StarAlgebras.AlgebraElement(a, RG) isa StarAlgebras.AlgebraElement
    @test all(RG(g) isa StarAlgebras.AlgebraElement{typeof(RG)} for g in b)

    @test_throws AssertionError StarAlgebras.AlgebraElement([1,2,3], RG)
    @test StarAlgebras.AlgebraElement([1,2,3,0,0,0,0,0,0,0,0,0,0], RG) isa StarAlgebras.AlgebraElement

    p = Word(A, [2,3])
    a = RG(p)
    @test StarAlgebras.coeffs(a)[b[p]] == 1
    @test StarAlgebras.coeffs(a) isa SparseVector
    @test all(StarAlgebras.coeffs(a)[i] == 0 for i in 1:length(b) if i ≠ b[p])
    @test a(p) == 1
    @test all(a(g) == 0 for g in b if g != p)

    @test sprint(show, zero(RG)) == "0·(id)"
    @test sprint(show, one(RG)) == "1·(id)"
    @test isone(one(a))
    @test iszero(zero(a))
    @test sprint(show, a) == "1·b·c"
    @test sprint(show, -a) == "-1·b·c"

    @test hash(a) == hash(one(RG) + RG(p) - one(RG))

    z = zeros(l)
    z[b[p]] = 1
    @test StarAlgebras.AlgebraElement(z, RG) == a

    @test StarAlgebras.supp(a) == [p]
    @test StarAlgebras.supp_ind(a) == [b[p]]

    s = one(first(b))
    @test a(s) == 0

    a[s] = 2

    @test StarAlgebras.coeffs(a)[b[s]] == 2
    @test a[b[s]] == 2
    @test a(s) == 2

    dense_a = StarAlgebras.AlgebraElement(Vector(StarAlgebras.coeffs(a)), RG)
    @test a == dense_a
    @test hash(a) == hash(dense_a)

    @test StarAlgebras.supp(a) == [s, p]
    @test StarAlgebras.supp_ind(a) == [b[s], b[p]] == StarAlgebras.supp_ind(dense_a)

    @test sprint(show, a) == "2·(id) +1·b·c"
    @test sprint(show, -a) == "-2·(id) -1·b·c"
    z[b[s]] = 2
    @test StarAlgebras.AlgebraElement(z, RG) == a
    @test sprint(show, StarAlgebras.AlgebraElement(z, RG)) == "2.0·(id) +1.0·b·c"

    @test LinearAlgebra.norm(a, 1) == 3
end