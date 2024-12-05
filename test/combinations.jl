# This section tests:
# - All supported combinations of integral(f, ::Geometry, ::IntegrationAlgorithm)
# - Invalid applications of integral aliases produce a descriptive error

#===============================================================================
                        Test Generation Infrastructure
===============================================================================#

@testsnippet Combinations begin
    using LinearAlgebra: norm
    using Meshes
    using MeshIntegrals
    using Unitful

    struct Callable{F <: Function}
        f::F
    end
    (c::Callable)(p) = c.f(p)

    # Stores a testable combination
    struct TestableGeometry{F <: Function, G <: Geometry, U <: Unitful.Quantity}
        integrand::F
        geometry::G
        solution::U
    end

    # Indicates which functions/rules are supported for a particular geometry
    struct SupportStatus
        lineintegral::Bool
        surfaceintegral::Bool
        volumeintegral::Bool
        gausskronrod::Bool
        gausslegendre::Bool
        hadaptivecubature::Bool
    end

    # Shortcut constructor for geometries with typical support structure
    function SupportStatus(sym::Symbol)
        if sym == :line
            aliases = Bool.((1, 0, 0))
            rules = Bool.((1, 1, 1))
            return SupportStatus(aliases..., rules...)
        elseif sym == :surface
            aliases = Bool.((0, 1, 0))
            rules = Bool.((1, 1, 1))
            return SupportStatus(aliases..., rules...)
        elseif sym == :volume
            aliases = Bool.((0, 0, 1))
            rules = Bool.((0, 1, 1))
            return SupportStatus(aliases..., rules...)
        else
            error("Unrecognized SupportStatus shortcut $(string(sym))")
        end
    end

    function runtests(testable::TestableGeometry, supports::SupportStatus)
        # Test alias functions
        for alias in (lineintegral, surfaceintegral, volumeintegral)
            alias_symbol = first(methods(alias)).name
            if getfield(supports, alias_symbol)
                @test alias(testable.integrand, testable.geometry) ≈ testable.solution
            else
                @test_throws "not supported" alias(testable.integrand, testable.geometry)
            end
        end

        iter_rules = (
            (supports.gausskronrod, GaussKronrod()),
            (supports.gausslegendre, GaussLegendre(100)),
            (supports.hadaptivecubature, HAdaptiveCubature())
        )

        # Test rules
        for (supported, rule) in iter_rules
            if supported
                # Scalar integrand
                sol = testable.solution
                @test integral(testable.integrand, testable.geometry, rule) ≈ sol

                # Callable integrand
                f = Callable(testable.integrand)
                @test integral(f, testable.geometry, rule) ≈ sol

                # Vector integrand
                fv(p) = fill(testable.integrand(p), 3)
                sol_v = fill(testable.solution, 3)
                @test integral(fv, testable.geometry, rule) ≈ sol_v
            else
                @test_throws "not supported" integral(testable.integrand, testable.geometry, rule)
            end # if
        end # for
    end # function
end #testsnippet

#===============================================================================
                         Create and Test Geometries
===============================================================================#

@testitem "Meshes.Ball 2D" setup=[Combinations] begin
    # Geometry
    origin = Point(0, 0)
    radius = 2.8
    ball = Ball(origin, radius)

    # Integrand & Solution
    function integrand(p::P) where {P <: Meshes.Point}
        r = ustrip(u"m", norm(to(p)))
        exp(-r^2)
    end
    solution = (π - π * exp(-radius^2)) * u"m^2"

    # Package and run tests
    testable = TestableGeometry(integrand, ball, solution)
    runtests(testable, SupportStatus(:surface))
end

@testitem "Meshes.Ball 3D" setup=[Combinations] begin
    using SpecialFunctions: erf

    # Geometry
    center = Point(1, 2, 3)
    radius = 2.8u"m"
    ball = Ball(center, radius)

    # Integrand & Solution
    function integrand(p::P) where {P <: Meshes.Point}
        offset = p - center
        ur = norm(offset)
        r = ustrip(u"m", ur)
        exp(-r^2)
    end
    solution = let r = ustrip(u"m", radius)
        (π^(3 / 2) * erf(r) - 2π * exp(-r^2) * r) * u"m^3"
    end

    # Package and run tests
    testable = TestableGeometry(integrand, ball, solution)
    runtests(testable, SupportStatus(:volume))
end

@testitem "Meshes.BezierCurve" setup=[Setup] begin
    # Geometry
    curve = BezierCurve([Point(t, sin(t), 0) for t in range(-π, π, length = 361)])

    # Integrand
    function f(p::P) where {P <: Meshes.Point}
        ux = ustrip(p.coords.x)
        (1 / sqrt(1 + cos(ux)^2)) * u"Ω/m"
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 2π * u"Ω"
    @test integral(f, curve, GaussLegendre(100))≈sol rtol=0.5e-2
    @test integral(f, curve, GaussKronrod())≈sol rtol=0.5e-2
    @test integral(f, curve, HAdaptiveCubature())≈sol rtol=0.5e-2

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, curve, GaussLegendre(100))≈vsol rtol=0.5e-2
    @test integral(fv, curve, GaussKronrod())≈vsol rtol=0.5e-2
    @test integral(fv, curve, HAdaptiveCubature())≈vsol rtol=0.5e-2

    # Integral aliases
    @test lineintegral(f, curve)≈sol rtol=0.5e-2
    @test_throws "not supported" surfaceintegral(f, curve)
    @test_throws "not supported" volumeintegral(f, curve)

    # Check Bezier-specific jacobian bounds
    @test_throws DomainError jacobian(curve, [1.1])
end

@testitem "Meshes.Box 1D" setup=[Setup] begin
    a = π
    box = Box(Point(0), Point(a))

    struct Integrand end
    fc = Integrand()

    function (::Integrand)(p::Meshes.Point)
        t = ustrip(p.coords.x)
        sqrt(a^2 - t^2) * u"Ω/m"
    end
    fcv(p) = fill(fc(p), 3)

    # Scalar integrand
    sol = π * a^2 / 4 * u"Ω"
    @test integral(fc, box, GaussLegendre(100))≈sol rtol=1e-6
    @test integral(fc, box, GaussKronrod()) ≈ sol
    @test integral(fc, box, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fcv, box, GaussLegendre(100))≈vsol rtol=1e-6
    @test integral(fcv, box, GaussKronrod()) ≈ vsol
    @test integral(fcv, box, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test lineintegral(fc, box) ≈ sol
    @test_throws "not supported" surfaceintegral(fc, box)
    @test_throws "not supported" volumeintegral(fc, box)
end

@testitem "Meshes.Box 2D" setup=[Setup] begin
    a = π
    box = Box(Point(0, 0), Point(a, a))

    function f(p::P) where {P <: Meshes.Point}
        x, y = ustrip.((p.coords.x, p.coords.y))
        (sqrt(a^2 - x^2) + sqrt(a^2 - y^2)) * u"Ω/m^2"
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 2a * (π * a^2 / 4) * u"Ω"
    @test integral(f, box, GaussLegendre(100))≈sol rtol=1e-6
    @test integral(f, box, GaussKronrod()) ≈ sol
    @test integral(f, box, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, box, GaussLegendre(100))≈vsol rtol=1e-6
    @test integral(fv, box, GaussKronrod()) ≈ vsol
    @test integral(fv, box, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, box)
    @test surfaceintegral(f, box) ≈ sol
    @test_throws "not supported" volumeintegral(f, box)

    # Test jacobian with wrong number of parametric coordinates
    @test_throws ArgumentError jacobian(box, zeros(3))
end

@testitem "Meshes.Box 3D" setup=[Setup] begin
    a = π
    box = Box(Point(0, 0, 0), Point(a, a, a))

    function f(p::P) where {P <: Meshes.Point}
        x, y, z = ustrip.((p.coords.x, p.coords.y, p.coords.z))
        (sqrt(a^2 - x^2) + sqrt(a^2 - y^2) + sqrt(a^2 - z^2)) * u"Ω/m^3"
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 3a^2 * (π * a^2 / 4) * u"Ω"
    @test integral(f, box, GaussLegendre(100))≈sol rtol=1e-6
    @test_throws "not supported" integral(f, box, GaussKronrod())
    @test integral(f, box, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, box, GaussLegendre(100))≈vsol rtol=1e-6
    @test_throws "not supported" integral(fv, box, GaussKronrod())
    @test integral(fv, box, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, box)
    @test_throws "not supported" surfaceintegral(f, box)
    @test volumeintegral(f, box) ≈ sol
end

@testitem "Meshes.Box 4D" tags=[:extended] setup=[Setup] begin
    a = π
    box = Box(Point(0, 0, 0, 0), Point(a, a, a, a))

    function f(p::P) where {P <: Meshes.Point}
        x1, x2, x3, x4 = ustrip.(to(p).coords)
        σ(x) = sqrt(a^2 - x^2)
        (σ(x1) + σ(x2) + σ(x3) + σ(x4)) * u"Ω/m^4"
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 4a^3 * (π * a^2 / 4) * u"Ω"
    @test integral(f, box, GaussLegendre(100))≈sol rtol=1e-6
    @test_throws "not supported" integral(f, box, GaussKronrod())
    @test integral(f, box, HAdaptiveCubature(rtol = 1e-6))≈sol rtol=1e-6

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, box, GaussLegendre(100))≈vsol rtol=1e-6
    @test_throws "not supported" integral(fv, box, GaussKronrod())
    @test integral(fv, box, HAdaptiveCubature(rtol = 1e-6))≈vsol rtol=1e-6

    # Integral aliases
    @test_throws "not supported" lineintegral(f, box)
    @test_throws "not supported" surfaceintegral(f, box)
    @test_throws "not supported" volumeintegral(f, box)
end

@testitem "Meshes.Circle" setup=[Combinations] begin
    # Geometry
    center = Point(1, 2, 3)
    n̂ = Vec(1 / 2, 1 / 2, sqrt(2) / 2)
    plane = Plane(center, n̂)
    radius = 4.4
    circle = Circle(plane, radius)

    # Integrand
    function integrand(p::P) where {P <: Meshes.Point}
        offset = p - center
        r = ustrip(u"m", norm(offset))
        exp(-r^2)
    end

    # Scalar integrand
    solution = 2π * radius * exp(-radius^2) * u"m"

    # Package and run tests
    testable = TestableGeometry(integrand, circle, solution)
    runtests(testable, SupportStatus(:line))
end

@testitem "Meshes.Cone" setup=[Combinations] begin
    # Geometry
    r = 2.5u"m"
    h = 3.5u"m"
    origin = Point(0, 0, 0)
    xy_plane = Plane(origin, Vec(0, 0, 1))
    base = Disk(xy_plane, r)
    apex = Point(0.0u"m", 0.0u"m", h)
    cone = Cone(base, apex)

    # Integrand
    integrand(p) = 1.0u"A"

    # Solution
    solution = (π * r^2 * h / 3) * u"A"

    # Package and run tests
    testable = TestableGeometry(integrand, cone, solution)
    runtests(testable, SupportStatus(:volume))
end

@testitem "Meshes.ConeSurface" setup=[Setup] begin
    r = 2.5u"m"
    h = 2.5u"m"
    origin = Point(0, 0, 0)
    xy_plane = Plane(origin, Vec(0, 0, 1))
    base = Disk(xy_plane, r)
    apex = Point(0.0u"m", 0.0u"m", h)
    cone = ConeSurface(base, apex)

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    _area_cone_rightcircular(h, r) = (π * r^2) + (π * r * hypot(h, r))

    # Scalar integrand
    sol = _area_cone_rightcircular(h, r)
    @test integral(f, cone, GaussLegendre(100))≈sol rtol=1e-6
    @test integral(f, cone, GaussKronrod())≈sol rtol=1e-6
    @test integral(f, cone, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, cone, GaussLegendre(100))≈vsol rtol=1e-6
    @test integral(fv, cone, GaussKronrod())≈vsol rtol=1e-6
    @test integral(fv, cone, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, cone)
    @test surfaceintegral(f, cone) ≈ sol
    @test_throws "not supported" volumeintegral(f, cone)
end

@testitem "Meshes.Cylinder" setup=[Combinations] begin
    # Geometry
    pt_w = Point(-1, 0, 0)
    pt_e = Point(1, 0, 0)
    cyl = Cylinder(pt_e, pt_w, 2.5)

    # Integrand & Solution
    integrand(p) = 1.0u"A"
    solution = Meshes.measure(cyl) * u"A"

    # Package and run tests
    testable = TestableGeometry(integrand, cyl, solution)
    runtests(testable, SupportStatus(:volume))
end

@testitem "Meshes.CylinderSurface" setup=[Combinations] begin
    # Geometry
    pt_w = Point(-1, 0, 0)
    pt_e = Point(1, 0, 0)
    cyl = CylinderSurface(pt_e, pt_w, 2.5)

    # Integrand & Solution
    integrand(p) = 1.0u"A"
    solution = Meshes.measure(cyl) * u"A"

    # Package and run tests
    testable = TestableGeometry(integrand, cyl, solution)
    runtests(testable, SupportStatus(:surface))
end

@testitem "Meshes.Disk" setup=[Combinations] begin
    # Geometry
    center = Point(1, 2, 3)
    n̂ = Vec(1 / 2, 1 / 2, sqrt(2) / 2)
    plane = Plane(center, n̂)
    radius = 2.5
    disk = Disk(plane, radius)

    # Integrand & Solution
    function integrand(p::P) where {P <: Meshes.Point}
        offset = p - center
        r = ustrip(u"m", norm(offset))
        exp(-r^2) * u"A"
    end
    solution = (π - π * exp(-radius^2)) * u"A*m^2"

    # Package and run tests
    testable = TestableGeometry(integrand, disk, solution)
    runtests(testable, SupportStatus(:surface))
end

@testitem "Meshes.Ellipsoid" setup=[Setup] begin
    origin = Point(0, 0, 0)
    radii = (1.0, 2.0, 0.5)
    ellipsoid = Ellipsoid(radii, origin)

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    # Tolerances are higher due to `measure` being only an approximation
    # Scalar integrand
    sol = Meshes.measure(ellipsoid)
    @test integral(f, ellipsoid, GaussLegendre(100))≈sol rtol=1e-2
    @test integral(f, ellipsoid, GaussKronrod())≈sol rtol=1e-2
    @test integral(f, ellipsoid, HAdaptiveCubature())≈sol rtol=1e-2

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, ellipsoid, GaussLegendre(100))≈vsol rtol=1e-2
    @test integral(fv, ellipsoid, GaussKronrod())≈vsol rtol=1e-2
    @test integral(fv, ellipsoid, HAdaptiveCubature())≈vsol rtol=1e-2

    # Integral aliases
    @test_throws "not supported" lineintegral(f, ellipsoid)
    @test surfaceintegral(f, ellipsoid)≈sol rtol=1e-2
    @test_throws "not supported" volumeintegral(f, ellipsoid)
end

@testitem "Meshes.FrustumSurface" setup=[Setup] begin
    # Create a frustum whose radius halves at the top,
    # i.e. the bottom half of a cone by height
    r_bot = 2.5u"m"
    r_top = 1.25u"m"
    cone_h = 2π * u"m"
    origin = Point(0, 0, 0)
    z = Vec(0, 0, 1)
    plane_bot = Plane(origin, z)
    disk_bot = Disk(plane_bot, r_bot)
    center_top = Point(0.0u"m", 0.0u"m", 0.5cone_h)
    plane_top = Plane(center_top, z)
    disk_top = Disk(plane_top, r_top)
    frustum = FrustumSurface(disk_bot, disk_top)

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    _area_base(r) = π * r^2
    _area_cone_walls(h, r) = π * r * hypot(h, r)

    # Scalar integrand
    sol = let
        area_walls_projected = _area_cone_walls(cone_h, r_bot)
        area_walls_missing = _area_cone_walls(0.5cone_h, r_top)
        area_walls = area_walls_projected - area_walls_missing
        area_walls + _area_base(r_top) + _area_base(r_bot)
    end
    @test integral(f, frustum, GaussLegendre(100))≈sol rtol=1e-6
    @test integral(f, frustum, GaussKronrod())≈sol rtol=1e-6
    @test integral(f, frustum, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, frustum, GaussLegendre(100))≈vsol rtol=1e-6
    @test integral(fv, frustum, GaussKronrod())≈vsol rtol=1e-6
    @test integral(fv, frustum, HAdaptiveCubature()) ≈ vsol
end

@testitem "Meshes.Hexahedron" setup=[Setup] begin
    hexahedron = Hexahedron(Point(0, 0, 0), Point(2, 0, 0), Point(2, 2, 0),
        Point(0, 2, 0), Point(0, 0, 2), Point(1, 0, 2), Point(1, 1, 2), Point(0, 1, 2))

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = Meshes.measure(hexahedron)
    @test integral(f, hexahedron, GaussLegendre(100)) ≈ sol
    @test_throws "not supported" integral(f, hexahedron, GaussKronrod())≈sol
    @test integral(f, hexahedron, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, hexahedron, GaussLegendre(100)) ≈ vsol
    @test_throws "not supported" integral(fv, hexahedron, GaussKronrod())≈vsol
    @test integral(fv, hexahedron, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, hexahedron)
    @test_throws "not supported" surfaceintegral(f, hexahedron)
    @test volumeintegral(f, hexahedron) ≈ sol
end

@testitem "Meshes.Line" setup=[Combinations] begin
    # Geometry
    a = Point(0, 0, 0)
    b = Point(1, 1, 1)
    line = Line(a, b)

    # Integrand & solution
    function integrand(p::P) where {P <: Meshes.Point}
        r = ustrip(u"m", norm(to(p)))
        exp(-r^2) * u"A"
    end
    solution = sqrt(π) * u"A*m"

    # Package and run tests
    testable = TestableGeometry(integrand, line, solution)
    runtests(testable, SupportStatus(:line))
end

@testitem "Meshes.ParaboloidSurface" setup=[Combinations] begin
    origin = Point(0, 0, 0)
    parab = ParaboloidSurface(origin, 2.5, 4.15)

    # Integrand & Solution
    integrand(p) = 1.0u"A"
    solution = Meshes.measure(parab) * u"A"

    # Package and run tests
    testable = TestableGeometry(integrand, parab, solution)
    runtests(testable, SupportStatus(:surface))
end

@testitem "ParametrizedCurve" setup=[Combinations] begin
    # ParametrizedCurve has been added in Meshes v0.51.20
    # If the version is specified as minimal compat bound in the Project.toml, the downgrade test fails
    if pkgversion(Meshes) >= v"0.51.20"
        using CoordRefSystems: Polar

        # Geometries
        # Parametrize a circle centered on origin with specified radius
        radius = 4.4
        curve_cart = ParametrizedCurve(
            t -> Point(radius * cos(t), radius * sin(t)), (0.0, 2π))
        curve_polar = ParametrizedCurve(t -> Point(Polar(radius, t)), (0.0, 2π))

        # Integrand & Solution
        function integrand(p::P) where {P <: Meshes.Point}
            ur = norm(to(p))
            r = ustrip(u"m", ur)
            exp(-r^2) * u"A"
        end
        solution = 2π * radius * exp(-radius^2) * u"A*m"

        # Package and run tests
        testable_cart = TestableGeometry(integrand, curve_cart, solution)
        runtests(testable_cart, SupportStatus(:line))
        testable_polar = TestableGeometry(integrand, curve_cart, solution)
        runtests(testable_polar, SupportStatus(:line))
    end
end

@testitem "Meshes.Plane" setup=[Combinations] begin
    # Geometry
    p = Point(0.0u"m", 0.0u"m", 0.0u"m")
    v = Vec(0.0u"m", 0.0u"m", 1.0u"m")
    plane = Plane(p, v)

    # Integrand & Solution
    function integrand(p::P) where {P <: Meshes.Point}
        r = ustrip(u"m", norm(to(p)))
        exp(-r^2) * u"A"
    end
    solution = π * u"A*m^2"

    # Package and run tests
    testable = TestableGeometry(integrand, plane, solution)
    runtests(testable, SupportStatus(:surface))
end

@testitem "Meshes.Quadrangle" setup=[Combinations] begin
    using SpecialFunctions: erf

    # Geometry
    quadrangle = Quadrangle((-1.0, 0.0), (-1.0, 1.0), (1.0, 1.0), (1.0, 0.0))

    # Integrand & Solution
    function integrand(p::P) where {P <: Meshes.Point}
        r = ustrip(u"m", norm(to(p)))
        exp(-r^2) * u"A"
    end
    solution = 0.5 * π * erf(1)^2 * u"A*m^2"

    # Package and run tests
    testable = TestableGeometry(integrand, quadrangle, solution)
    runtests(testable, SupportStatus(:surface))
end

@testitem "Meshes.Ray" setup=[Combinations] begin
    # Geometry
    a = Point(0, 0, 0)
    v = Vec(1, 1, 1)
    ray = Ray(a, v)

    # Integrand & Solution
    function integrand(p::P) where {P <: Meshes.Point}
        r = ustrip(u"m", norm(to(p)))
        exp(-r^2) * u"A"
    end
    solution = sqrt(π) / 2 * u"A*m"

    # Package and run tests
    testable = TestableGeometry(integrand, ray, solution)
    runtests(testable, SupportStatus(:line))
end

@testitem "Meshes.Ring" setup=[Setup] begin
    pt_a = Point(0.0u"m", 0.0u"m", 0.0u"m")
    pt_b = Point(1.0u"m", 0.0u"m", 0.0u"m")
    pt_c = Point(1.0u"m", 1.0u"m", 0.0u"m")
    pt_d = Point(1.0u"m", 1.0u"m", 1.0u"m")
    rope = Ring(pt_a, pt_b, pt_c, pt_d, pt_c, pt_b)

    function f(p::P) where {P <: Meshes.Point}
        x, y, z = (p.coords.x, p.coords.y, p.coords.z)
        (x + 2y + 3z) * u"A/m^2"
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 14.0u"A"
    @test integral(f, rope, GaussLegendre(100)) ≈ sol
    @test integral(f, rope, GaussKronrod()) ≈ sol
    @test integral(f, rope, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, rope, GaussLegendre(100)) ≈ vsol
    @test integral(fv, rope, GaussKronrod()) ≈ vsol
    @test integral(fv, rope, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test lineintegral(f, rope) ≈ sol
    @test_throws "not supported" surfaceintegral(f, rope)
    @test_throws "not supported" volumeintegral(f, rope)
end

@testitem "Meshes.Rope" setup=[Setup] begin
    pt_a = Point(0.0u"m", 0.0u"m", 0.0u"m")
    pt_b = Point(1.0u"m", 0.0u"m", 0.0u"m")
    pt_c = Point(1.0u"m", 1.0u"m", 0.0u"m")
    pt_d = Point(1.0u"m", 1.0u"m", 1.0u"m")
    rope = Rope(pt_a, pt_b, pt_c, pt_d)

    function f(p::P) where {P <: Meshes.Point}
        x, y, z = (p.coords.x, p.coords.y, p.coords.z)
        (x + 2y + 3z) * u"A/m^2"
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 7.0u"A"
    @test integral(f, rope, GaussLegendre(100)) ≈ sol
    @test integral(f, rope, GaussKronrod()) ≈ sol
    @test integral(f, rope, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, rope, GaussLegendre(100)) ≈ vsol
    @test integral(fv, rope, GaussKronrod()) ≈ vsol
    @test integral(fv, rope, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test lineintegral(f, rope) ≈ sol
    @test_throws "not supported" surfaceintegral(f, rope)
    @test_throws "not supported" volumeintegral(f, rope)
end

@testitem "Meshes.Segment" setup=[Setup] begin
    # Connect a line segment from the origin to an arbitrary point on the unit sphere
    φ, θ = (7pi / 6, pi / 3)  # Arbitrary spherical angles
    pt_a = Point(0.0u"m", 0.0u"m", 0.0u"m")
    pt_b = Point(sin(θ) * cos(φ) * u"m", sin(θ) * sin(φ) * u"m", cos(θ) * u"m")
    segment = Segment(pt_a, pt_b)

    a, b = (7.1, 4.6)  # arbitrary constants > 0

    function f(p::P) where {P <: Meshes.Point}
        r = ustrip(u"m", norm(to(p)))
        exp(r * log(a) + (1 - r) * log(b))
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = (a - b) / (log(a) - log(b)) * u"m"
    @test integral(f, segment, GaussLegendre(100)) ≈ sol
    @test integral(f, segment, GaussKronrod()) ≈ sol
    @test integral(f, segment, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, segment, GaussLegendre(100)) ≈ vsol
    @test integral(fv, segment, GaussKronrod()) ≈ vsol
    @test integral(fv, segment, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test lineintegral(f, segment) ≈ sol
    @test_throws "not supported" surfaceintegral(f, segment)
    @test_throws "not supported" volumeintegral(f, segment)
end

@testitem "Meshes.Sphere 2D" setup=[Setup] begin
    origin = Point(0, 0)
    radius = 4.4
    sphere = Sphere(origin, radius)

    function f(p::P) where {P <: Meshes.Point}
        r = ustrip(u"m", norm(to(p)))
        exp(-r^2)
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 2π * radius * exp(-radius^2) * u"m"
    @test integral(f, sphere, GaussLegendre(100)) ≈ sol
    @test integral(f, sphere, GaussKronrod()) ≈ sol
    @test integral(f, sphere, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, sphere, GaussLegendre(100)) ≈ vsol
    @test integral(fv, sphere, GaussKronrod()) ≈ vsol
    @test integral(fv, sphere, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test lineintegral(f, sphere) ≈ sol
    @test_throws "not supported" surfaceintegral(f, sphere)
    @test_throws "not supported" volumeintegral(f, sphere)
end

@testitem "Meshes.Sphere 3D" setup=[Setup] begin
    using CoordRefSystems: Cartesian, Spherical

    center = Point(1, 2, 3)
    radius = 4.4u"m"
    sphere = Sphere(center, radius)

    function f(p::P) where {P <: Meshes.Point}
        rθφ = convert(Spherical, Cartesian((p - center)...))
        r = ustrip(rθφ.r)
        θ = ustrip(rθφ.θ)
        φ = ustrip(rθφ.ϕ)
        sin(φ)^2 + cos(θ)^2
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = (2π * radius^2) + (4π / 3 * radius^2)
    @test integral(f, sphere, GaussLegendre(100)) ≈ sol
    @test integral(f, sphere, GaussKronrod()) ≈ sol
    @test integral(f, sphere, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, sphere, GaussLegendre(100)) ≈ vsol
    @test integral(fv, sphere, GaussKronrod()) ≈ vsol
    @test integral(fv, sphere, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, sphere)
    @test surfaceintegral(f, sphere) ≈ sol
    @test_throws "not supported" volumeintegral(f, sphere)
end

@testitem "Meshes.Tetrahedron" setup=[Setup] begin
    pt_n = Point(0, 3, 0)
    pt_w = Point(-7, 0, 0)
    pt_e = Point(8, 0, 0)
    ẑ = Vec(0, 0, 1)
    tetrahedron = Tetrahedron(pt_n, pt_w, pt_e, pt_n + ẑ)

    f(p) = 1.0u"A"
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = Meshes.measure(tetrahedron) * u"A"
    @test integral(f, tetrahedron, GaussLegendre(100)) ≈ sol
    @test_throws "not supported" integral(f, tetrahedron, GaussKronrod())
    @test integral(f, tetrahedron, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, tetrahedron, GaussLegendre(100)) ≈ vsol
    @test_throws "not supported" integral(fv, tetrahedron, GaussKronrod())
    @test integral(fv, tetrahedron, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, tetrahedron)
    @test_throws "not supported" surfaceintegral(f, tetrahedron)
    @test volumeintegral(f, tetrahedron, HAdaptiveCubature()) ≈ sol
end

@testitem "Meshes.Torus" setup=[Setup] begin
    origin = Point(0, 0, 0)
    ẑ = Vec(0, 0, 1)
    torus = Torus(origin, ẑ, 3.5, 1.25)

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = Meshes.measure(torus)
    @test integral(f, torus, GaussLegendre(100)) ≈ sol
    @test integral(f, torus, GaussKronrod()) ≈ sol
    @test integral(f, torus, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, torus, GaussLegendre(100)) ≈ vsol
    @test integral(fv, torus, GaussKronrod()) ≈ vsol
    @test integral(fv, torus, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, torus)
    @test surfaceintegral(f, torus) ≈ sol
    @test_throws "not supported" volumeintegral(f, torus)
end

@testitem "Meshes.Triangle" setup=[Setup] begin
    pt_n = Point(0, 1, 0)
    pt_w = Point(-1, 0, 0)
    pt_e = Point(1, 0, 0)
    triangle = Triangle(pt_e, pt_n, pt_w)

    f(p) = 1.0u"A"
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = Meshes.measure(triangle) * u"A"
    @test integral(f, triangle, GaussLegendre(100)) ≈ sol
    @test integral(f, triangle, GaussKronrod()) ≈ sol
    @test integral(f, triangle, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, triangle, GaussLegendre(100)) ≈ vsol
    @test integral(fv, triangle, GaussKronrod()) ≈ vsol
    @test integral(fv, triangle, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, triangle)
    @test surfaceintegral(f, triangle) ≈ sol
    @test_throws "not supported" volumeintegral(f, triangle)
end
