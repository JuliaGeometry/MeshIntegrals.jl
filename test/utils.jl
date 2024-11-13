@testitem "Utilities" setup=[Setup] begin
    using LinearAlgebra: norm
    using MeshIntegrals: _units, _zeros, _ones

    # _KVector
    v = Meshes.Vec(3, 4)
    @test norm(MeshIntegrals._KVector(v)) ≈ 5.0u"m"

    # _units
    p = Point(1.0u"cm", 2.0u"mm", 3.0u"m")
    @test _units(p) == u"m"

    # _zeros
    @test _zeros(2) == (0.0, 0.0)
    @test _zeros(Float32, 2) == (0.0f0, 0.0f0)

    # _ones
    @test _ones(2) == (1.0, 1.0)
    @test _ones(Float32, 2) == (1.0f0, 1.0f0)
end

@testitem "DifferentiationMethod" setup=[Setup] begin
    using MeshIntegrals: _has_analytical, _default_method

    # _has_analytical of instances
    bezier = BezierCurve([Point(t, sin(t), 0.0) for t in range(-π, π, length = 361)])
    @test _has_analytical(bezier) == false
    sphere = Sphere(Point(0, 0, 0), 1.0)
    @test _has_analytical(sphere) == false

    # _has_analytical of types
    @test _has_analytical(Meshes.BezierCurve) == false
    @test _has_analytical(Meshes.Sphere) == false

    # _default_method
    @test _default_method(Meshes.BezierCurve) isa FiniteDifference
    @test _default_method(bezier) isa FiniteDifference
    @test _default_method(Meshes.Sphere) isa FiniteDifference
    @test _default_method(sphere) isa FiniteDifference

    # FiniteDifference
    @test FiniteDifference().ε ≈ 1e-6
end

@testitem "_ParametricGeometry" setup=[Setup] begin
    using MeshIntegrals: _ParametricGeometry, _parametric

    # paramdim(::_ParametricGeometry)
    segment = Segment(Point(0, 0), Point(1, 1))
    f(t) = segment(t)
    geometry = _ParametricGeometry(f, 1)
    @test paramdim(geometry) == 1

    # _parametric bounds checks
    triangle = Triangle(Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1))
    @test_throws DomainError _parametric(triangle, 1.1, 0.0)
    tetrahedron = let
        a = Point(0, 3, 0)
        b = Point(-7, 0, 0)
        c = Point(8, 0, 0)
        ẑ = Vec(0, 0, 1)
        Tetrahedron(a, b, c, a + ẑ)
    end
    @test_throws DomainError _parametric(tetrahedron, 1.1, 0.0, 0.0)
end
