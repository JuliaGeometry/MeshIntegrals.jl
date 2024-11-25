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
    using MeshIntegrals: _has_analytical, _default_method, _guarantee_analytical

    # _has_analytical of instances
    bezier = BezierCurve([Point(t, sin(t), 0.0) for t in range(-π, π, length = 361)])
    @test _has_analytical(bezier) == true
    sphere = Sphere(Point(0, 0, 0), 1.0)
    @test _has_analytical(sphere) == false

    # _has_analytical of types
    @test _has_analytical(Meshes.BezierCurve) == true
    @test _has_analytical(Meshes.Line) == true
    @test _has_analytical(Meshes.Plane) == true
    @test _has_analytical(Meshes.Ray) == true
    @test _has_analytical(Meshes.Sphere) == false
    @test _has_analytical(Meshes.Tetrahedron) == true
    @test _has_analytical(Meshes.Triangle) == true

    # _guarantee_analytical
    @test _guarantee_analytical(Meshes.Line, Analytical()) === nothing
    @test_throws "Analytical" _guarantee_analytical(Meshes.Line, FiniteDifference())

    # _default_method
    @test _default_method(Meshes.BezierCurve) isa Analytical
    @test _default_method(bezier) isa Analytical
    @test _default_method(Meshes.Sphere) isa FiniteDifference
    @test _default_method(sphere) isa FiniteDifference

    # FiniteDifference
    @test FiniteDifference().ε ≈ 1e-6
end
