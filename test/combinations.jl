"""
This file includes tests for:
- Supported combinations {::Geometry, ::IntegrationRule} for `integral(f, geometry, rule)`.
- Integrand functions returning a `Unitful.Quantity` or a `Vector{Unitful.Quantity}`.
- Callable objects in place of an integrand function.
- Unsupported combinations produce a descriptive and useful error message.
- The appropriate alias function, e.g. `lineintegral` for 1D geometries, works.
- Invalid applications of integral aliases produce a descriptive and useful error message.
- (Planned) Usage of non-default DifferentiationMethods.
"""

#===============================================================================
                        Test Generation Infrastructure
===============================================================================#

@testsnippet Combinations begin
    using CoordRefSystems
    import LinearAlgebra: norm
    using Meshes
    using MeshIntegrals
    import Unitful: @u_str, Quantity, ustrip
    import Enzyme

    # Used for testing callable objects as integrand functions
    struct Callable{F <: Function}
        f::F
    end
    (c::Callable)(p::Meshes.Point) = c.f(p)

    # Stores a testable combination
    struct TestableGeometry{F <: Function, G <: Meshes.GeometryOrDomain, U <: Quantity}
        integrand::F
        geometry::G
        solution::U
    end

    # Used to indicate which features are supported for a particular geometry
    struct SupportStatus
        # Alias Functions
        lineintegral::Bool
        surfaceintegral::Bool
        volumeintegral::Bool
        # IntegrationRules
        gausskronrod::Bool
        gausslegendre::Bool
        hadaptivecubature::Bool
        # DifferentiationMethods
        autoenzyme::Bool
    end

    # Shortcut constructor for geometries with typical support structure
    function SupportStatus(geometry::G) where {G <: Meshes.GeometryOrDomain}
        # Check whether AutoEnzyme should be supported, i.e. not on blacklist
        unsupported_Gs = Union{BezierCurve, Cylinder, CylinderSurface, ParametrizedCurve}
        autoenzyme = !(G <: unsupported_Gs)

        N = Meshes.paramdim(geometry)
        if N == 1
            # line/curve
            aliases = Bool.((1, 0, 0))
            rules = Bool.((1, 1, 1))
            return SupportStatus(aliases..., rules..., autoenzyme)
        elseif N == 2
            # surface
            aliases = Bool.((0, 1, 0))
            rules = Bool.((1, 1, 1))
            return SupportStatus(aliases..., rules..., autoenzyme)
        elseif N == 3
            # volume
            aliases = Bool.((0, 0, 1))
            rules = Bool.((0, 1, 1))
            return SupportStatus(aliases..., rules..., autoenzyme)
        else
            # ≥4D
            aliases = Bool.((0, 0, 0))
            rules = Bool.((0, 1, 1))
            return SupportStatus(aliases..., rules..., autoenzyme)
        end #if
    end # function

    # Generate applicable tests for this geometry
    function runtests(testable::TestableGeometry; rtol = sqrt(eps()))
        # Determine support matrix for this geometry
        supports = SupportStatus(testable.geometry)

        # Ensure consistency of SupportStatus with supports_autoenzyme
        @test MeshIntegrals.supports_autoenzyme(testable.geometry) == supports.autoenzyme

        # Test alias functions
        for alias in (lineintegral, surfaceintegral, volumeintegral)
            # if supports.alias
            if getfield(supports, nameof(alias))
                @test alias(testable.integrand, testable.geometry)≈testable.solution rtol=rtol
            else
                @test_throws "not supported" alias(testable.integrand, testable.geometry)
            end
        end # for

        # Iteratively test all IntegrationRules
        iter_rules = (
            (supports.gausskronrod, GaussKronrod()),
            (supports.gausslegendre, GaussLegendre(100)),
            (supports.hadaptivecubature, HAdaptiveCubature())
        )
        for (supported, rule) in iter_rules
            if supported
                # Scalar integrand
                sol = testable.solution
                @test integral(testable.integrand, testable.geometry, rule)≈sol rtol=rtol

                # Callable integrand
                f = Callable(testable.integrand)
                @test integral(f, testable.geometry, rule)≈sol rtol=rtol

                # Vector integrand
                fv(p) = fill(testable.integrand(p), 3)
                sol_v = fill(testable.solution, 3)
                @test integral(fv, testable.geometry, rule)≈sol_v rtol=rtol
            else
                f = testable.integrand
                geometry = testable.geometry
                @test_throws "not supported" integral(f, geometry, rule)
            end
        end # for

        # Iteratively test all DifferentiationMethods
        iter_diff_methods = (
            (true, FiniteDifference()),
            (supports.autoenzyme, AutoEnzyme())
        )
        for (supported, method) in iter_diff_methods
            # Aliases for improved code readability
            f = testable.integrand
            geometry = testable.geometry
            sol = testable.solution

            if supported
                @test integral(f, geometry; diff_method = method)≈sol rtol=rtol
            else
                @test_throws "not supported" integral(f, geometry; diff_method = method)
            end
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
    function integrand(p::Meshes.Point)
        r = ustrip(u"m", norm(to(p)))
        exp(-r^2) * u"A"
    end
    solution = (π - π * exp(-radius^2)) * u"A*m^2"

    # Package and run tests
    testable = TestableGeometry(integrand, ball, solution)
    runtests(testable)
end

@testitem "Meshes.Ball 3D" setup=[Combinations] begin
    using SpecialFunctions: erf

    # Geometry
    center = Point(1, 2, 3)
    radius = 2.8u"m"
    r = ustrip(u"m", radius)
    ball = Ball(center, radius)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        offset = p - center
        ur = norm(offset)
        r = ustrip(u"m", ur)
        exp(-r^2) * u"A"
    end
    solution = (π^(3 / 2) * erf(r) - 2π * exp(-r^2) * r) * u"A*m^3"

    # Package and run tests
    testable = TestableGeometry(integrand, ball, solution)
    runtests(testable)
end

@testitem "Meshes.BezierCurve" setup=[Combinations] begin
    # Geometry
    curve = BezierCurve([Point(t, sin(t), 0) for t in range(-π, π, length = 361)])

    # Integrand
    function integrand(p::Meshes.Point)
        ux = ustrip(p.coords.x)
        (1 / sqrt(1 + cos(ux)^2)) * u"Ω"
    end
    solution = 2π * u"Ω*m"

    # Package and run tests
    testable = TestableGeometry(integrand, curve, solution)
    runtests(testable; rtol = 0.5e-2)
end

@testitem "Meshes.Box 1D" setup=[Combinations] begin
    # Geometry
    a = π
    box = Box(Point(0), Point(a))

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        x₁ = only(ustrip.(to(p)))
        √(a^2 - x₁^2) * u"A"
    end
    solution = π * a^2 / 4 * u"A*m"

    # Package and run tests
    testable = TestableGeometry(integrand, box, solution)
    runtests(testable; rtol = 1e-6)
end

@testitem "Meshes.Box 2D" setup=[Combinations] begin
    # Geometry
    a = π
    box = Box(Point(0, 0), Point(a, a))

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        x₁, x₂ = ustrip.(to(p))
        (√(a^2 - x₁^2) + √(a^2 - x₂^2)) * u"A"
    end
    solution = 2a * (π * a^2 / 4) * u"A*m^2"

    # Package and run tests
    testable = TestableGeometry(integrand, box, solution)
    runtests(testable; rtol = 1e-6)
end

@testitem "Meshes.Box 3D" setup=[Combinations] begin
    # Geometry
    a = π
    box = Box(Point(0, 0, 0), Point(a, a, a))

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        x₁, x₂, x₃ = ustrip.(to(p))
        (√(a^2 - x₁^2) + √(a^2 - x₂^2) + √(a^2 - x₃^2)) * u"A"
    end
    solution = 3a^2 * (π * a^2 / 4) * u"A*m^3"

    # Package and run tests
    testable = TestableGeometry(integrand, box, solution)
    runtests(testable; rtol = 1e-6)
end

@testitem "Meshes.Box 4D" tags=[:extended] setup=[Combinations] begin
    # Geometry
    a = π
    box = Box(Point(0, 0, 0, 0), Point(a, a, a, a))

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        x₁, x₂, x₃, x₄ = ustrip.((to(p)))
        (√(a^2 - x₁^2) + √(a^2 - x₂^2) + √(a^2 - x₃^2) + √(a^2 - x₄^2)) * u"A"
    end
    solution = 4a^3 * (π * a^2 / 4) * u"A*m^4"

    # Package and run tests
    testable = TestableGeometry(integrand, box, solution)
    runtests(testable; rtol = 1e-6)
end

@testitem "Meshes.CartesianGrid" setup=[Combinations] begin
    # Geometry
    a = π
    start = Point(0, 0)
    finish = Point(a, a)
    dims = (4, 4)
    grid = CartesianGrid(start, finish, dims = dims)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        x₁, x₂ = ustrip.(to(p))
        (√(a^2 - x₁^2) + √(a^2 - x₂^2)) * u"A"
    end
    solution = 2a * (π * a^2 / 4) * u"A*m^2"

    # Package and run tests
    testable = TestableGeometry(integrand, grid, solution)
    runtests(testable; rtol = 1e-6)
end

@testitem "Meshes.Circle" setup=[Combinations] begin
    # Geometry
    center = Point(1, 2, 3)
    n̂ = Vec(1 / 2, 1 / 2, sqrt(2) / 2)
    plane = Plane(center, n̂)
    radius = 4.4
    circle = Circle(plane, radius)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        offset = p - center
        r = ustrip(u"m", norm(offset))
        exp(-r^2) * u"A"
    end
    solution = 2π * radius * exp(-radius^2) * u"A*m"

    # Package and run tests
    testable = TestableGeometry(integrand, circle, solution)
    runtests(testable)
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

    # Integrand & Solution
    integrand(p) = 1.0u"A"
    solution = (π * r^2 * h / 3) * u"A"

    # Package and run tests
    testable = TestableGeometry(integrand, cone, solution)
    runtests(testable)
end

@testitem "Meshes.ConeSurface" setup=[Combinations] begin
    # Geometry
    r = 2.5u"m"
    h = 3.5u"m"
    origin = Point(0, 0, 0)
    xy_plane = Plane(origin, Vec(0, 0, 1))
    base = Disk(xy_plane, r)
    apex = Point(0.0u"m", 0.0u"m", h)
    cone = ConeSurface(base, apex)

    # Integrand & Solution
    integrand(p) = 1.0u"A"
    solution = ((π * r^2) + (π * r * hypot(h, r))) * u"A"

    # Package and run tests
    testable = TestableGeometry(integrand, cone, solution)
    runtests(testable; rtol = 1e-6)
end

@testitem "Meshes.Cylinder" setup=[Combinations] begin
    # Geometry
    h = 8.5u"m"
    ρ₀ = 1.3u"m"
    pt_a = Point(0u"m", 0u"m", 0u"m")
    pt_b = Point(0u"m", 0u"m", h)
    cyl = Cylinder(pt_a, pt_b, ρ₀)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        p_cyl = convert(Cylindrical, Cartesian(to(p)...))
        ρ = p_cyl.ρ
        φ = p_cyl.ϕ
        z = p_cyl.z
        ρ^(-1) * (ρ + φ * u"m" + z) * u"A"
    end
    solution = ((π * h * ρ₀^2) + (π * h^2 * ρ₀) + (2π * π * u"m" * h * ρ₀)) * u"A"

    # Package and run tests
    testable = TestableGeometry(integrand, cyl, solution)
    runtests(testable)
end

@testitem "Meshes.CylinderSurface" setup=[Combinations] begin
    # Geometry
    h = 8.5u"m"
    ρ₀ = 1.3u"m"
    pt_a = Point(0u"m", 0u"m", 0u"m")
    pt_b = Point(0u"m", 0u"m", h)
    cyl = CylinderSurface(pt_a, pt_b, ρ₀)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        p_cyl = convert(Cylindrical, Cartesian(to(p)...))
        ρ = p_cyl.ρ
        φ = p_cyl.ϕ
        z = p_cyl.z
        ρ^(-1) * (ρ + φ * u"m" + z) * u"A"
    end
    solution = let
        disk_a = (2π * h * ρ₀) + (π * ρ₀^2) + (π * u"m" * ρ₀ * 2π)
        disk_b = (π * ρ₀^2) + (π * u"m" * ρ₀ * 2π)
        walls = (2π * h * ρ₀) + (2π^2 * u"m" * h) + (π * h^2)
        (disk_a + disk_b + walls) * u"A"
    end

    # Package and run tests
    testable = TestableGeometry(integrand, cyl, solution)
    runtests(testable)
end

@testitem "Meshes.Disk" setup=[Combinations] begin
    # Geometry
    center = Point(1, 2, 3)
    n̂ = Vec(1 / 2, 1 / 2, sqrt(2) / 2)
    plane = Plane(center, n̂)
    radius = 2.5
    disk = Disk(plane, radius)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        offset = p - center
        r = ustrip(u"m", norm(offset))
        exp(-r^2) * u"A"
    end
    solution = (π - π * exp(-radius^2)) * u"A*m^2"

    # Package and run tests
    testable = TestableGeometry(integrand, disk, solution)
    runtests(testable)
end

@testitem "Meshes.Ellipsoid" setup=[Combinations] begin
    # Geometry
    origin = Point(0, 0, 0)
    R = r₁ = r₂ = r₃ = 4.1
    ellipsoid = Ellipsoid((r₁, r₂, r₃), origin)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        x, y, z = ustrip.(u"m", to(p))
        (z^2) * u"A"
    end
    solution = (4π * R^4 / 3) * u"A*m^2"

    # Package and run tests
    # Tolerances are higher due to `measure` being only an approximation
    testable = TestableGeometry(integrand, ellipsoid, solution)
    runtests(testable; rtol = 1e-2)
end

@testitem "Meshes.Frustum" setup=[Combinations] begin
    if pkgversion(Meshes) >= v"0.52.12"
        # Geometry
        r = 2.5u"m"
        h = 3.5u"m"
        origin = Point(0, 0, 0)
        midpoint = Point(0.0u"m", 0.0u"m", 0.5 * h)
        apex = Point(0.0u"m", 0.0u"m", h)
        ẑ = Vec(0, 0, 1)
        xy_plane = Plane(origin, ẑ)
        base = Disk(xy_plane, r)
        mid_plane = Plane(midpoint, ẑ)
        mid_disk = Disk(mid_plane, 0.5 * r)
        frustum = Frustum(base, mid_disk)

        # Integrand & Solution
        integrand(p) = 1.0u"A"
        solution = (7 // 8) * (π * r^2 * h / 3) * u"A"

        # Package and run tests
        testable = TestableGeometry(integrand, frustum, solution)
        runtests(testable)
    end
end

@testitem "Meshes.FrustumSurface" setup=[Combinations] begin
    # Geometry
    # Create a frustum whose radius halves at the top, i.e. the bottom half of a cone
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

    # Integrand & Solution
    integrand(p) = 1.0u"A"
    _area_base(r) = π * r^2
    _area_cone_walls(h, r) = π * r * hypot(h, r)
    solution = let
        area_walls_projected = _area_cone_walls(cone_h, r_bot)
        area_walls_missing = _area_cone_walls(0.5cone_h, r_top)
        area_walls = area_walls_projected - area_walls_missing
        area_total = area_walls + _area_base(r_top) + _area_base(r_bot)
        area_total * u"A"
    end

    # Package and run tests
    testable = TestableGeometry(integrand, frustum, solution)
    runtests(testable; rtol = 1e-6)
end

@testitem "Meshes.Hexahedron" setup=[Combinations] begin
    # Geometry
    a = π
    box = Box(Point(0, 0, 0), Point(a, a, a))
    hexahedron = Hexahedron(discretize(box).vertices...)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        x₁, x₂, x₃ = ustrip.(to(p))
        (√(a^2 - x₁^2) + √(a^2 - x₂^2) + √(a^2 - x₃^2)) * u"A"
    end
    solution = 3a^2 * (π * a^2 / 4) * u"A*m^3"

    # Package and run tests
    testable = TestableGeometry(integrand, hexahedron, solution)
    runtests(testable; rtol = 1e-6)
end

@testitem "Meshes.Line" setup=[Combinations] begin
    # Geometry
    a = Point(0, 0, 0)
    b = Point(1, 1, 1)
    line = Line(a, b)

    # Integrand & solution
    function integrand(p::Meshes.Point)
        r = ustrip(u"m", norm(to(p)))
        exp(-r^2) * u"A"
    end
    solution = sqrt(π) * u"A*m"

    # Package and run tests
    testable = TestableGeometry(integrand, line, solution)
    runtests(testable)
end

@testitem "Meshes.ParaboloidSurface" setup=[Combinations] begin
    # Geometry
    origin = Point(0, 0, 0)
    parab = ParaboloidSurface(origin, 2.5, 4.15)

    # Integrand & Solution
    integrand(p) = 1.0u"A"
    solution = Meshes.measure(parab) * u"A"

    # Package and run tests
    testable = TestableGeometry(integrand, parab, solution)
    runtests(testable)
end

@testitem "Meshes.ParametrizedCurve" setup=[Combinations] begin
    using CoordRefSystems: Polar

    # Geometries
    # Parametrize a circle centered on origin with specified radius
    radius = 4.4
    curve_cart = ParametrizedCurve(
        t -> Point(radius * cos(t), radius * sin(t)), (0.0, 2π))
    curve_polar = ParametrizedCurve(t -> Point(Polar(radius, t)), (0.0, 2π))

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        ur = norm(to(p))
        r = ustrip(u"m", ur)
        exp(-r^2) * u"A"
    end
    solution = 2π * radius * exp(-radius^2) * u"A*m"

    # Package and run tests
    testable_cart = TestableGeometry(integrand, curve_cart, solution)
    runtests(testable_cart)
    testable_polar = TestableGeometry(integrand, curve_polar, solution)
    runtests(testable_polar)
end

@testitem "Meshes.Plane" setup=[Combinations] begin
    # Geometry
    p = Point(0.0u"m", 0.0u"m", 0.0u"m")
    v = Vec(0.0u"m", 0.0u"m", 1.0u"m")
    plane = Plane(p, v)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        r = ustrip(u"m", norm(to(p)))
        exp(-r^2) * u"A"
    end
    solution = π * u"A*m^2"

    # Package and run tests
    testable = TestableGeometry(integrand, plane, solution)
    runtests(testable)
end

@testitem "Meshes.PolyArea" setup=[Combinations] begin
    # Geometry
    a, b, c = 0.4, 0.6, 1
    outer = [(0, 0), (c, 0), (c, c), (0, c)]
    hole = [(a, a), (b, a), (b, b), (a, b)]
    area = PolyArea([outer, hole])

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        x, y = ustrip.(u"m", to(p))
        2x * u"A"
    end
    solution = (c^2 - (b - a) * (b^2 - a^2)) * u"A*m^2"

    # Package and run tests
    testable = TestableGeometry(integrand, area, solution)
    runtests(testable)
end

@testitem "Meshes.Pyramid" setup=[Combinations] begin
    if pkgversion(Meshes) >= v"0.52.12"
        # Geometry
        a = Point(-1, -1, 0)
        b = Point(1, -1, 0)
        c = Point(1, 1, 0)
        d = Point(-1, 1, 0)
        apex = Point(0, 0, 1)
        pyramid = Pyramid(a, b, c, d, apex)
        # Integrand & Solution
        integrand(p) = 1.0u"A"
        w = norm(b - a)
        h = norm(d - a)
        solution = (1 // 3) * w * h * u"A*m"
        # Package and run tests
        testable = TestableGeometry(integrand, pyramid, solution)
        runtests(testable)
    end
end

@testitem "Meshes.Quadrangle" setup=[Combinations] begin
    using SpecialFunctions: erf

    # Geometry
    quadrangle = Quadrangle((-1.0, 0.0), (-1.0, 1.0), (1.0, 1.0), (1.0, 0.0))

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        r = ustrip(u"m", norm(to(p)))
        exp(-r^2) * u"A"
    end
    solution = 0.5 * π * erf(1)^2 * u"A*m^2"

    # Package and run tests
    testable = TestableGeometry(integrand, quadrangle, solution)
    runtests(testable)
end

@testitem "Meshes.Ray" setup=[Combinations] begin
    # Geometry
    a = Point(0, 0, 0)
    v = Vec(1, 1, 1)
    ray = Ray(a, v)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        r = ustrip(u"m", norm(to(p)))
        exp(-r^2) * u"A"
    end
    solution = sqrt(π) / 2 * u"A*m"

    # Package and run tests
    testable = TestableGeometry(integrand, ray, solution)
    runtests(testable)
end

@testitem "Meshes.RegularGrid" setup=[Combinations] begin
    # Geometry
    a = π
    start = Point(0, 0)
    finish = Point(a, a)
    dims = (4, 4)
    grid = RegularGrid(start, finish, dims = dims)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        x₁, x₂ = ustrip.(to(p))
        (√(a^2 - x₁^2) + √(a^2 - x₂^2)) * u"A"
    end
    solution = 2a * (π * a^2 / 4) * u"A*m^2"

    # Package and run tests
    testable = TestableGeometry(integrand, grid, solution)
    runtests(testable; rtol = 1e-6)
end

@testitem "Meshes.Ring" setup=[Combinations] begin
    # Geometry
    a = Point(0, 0, 0)
    b = Point(1, 0, 0)
    c = Point(1, 1, 0)
    d = Point(1, 1, 1)
    ring = Ring(a, b, c, d, c, b)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        x, y, z = ustrip.((p.coords.x, p.coords.y, p.coords.z))
        (x + 2y + 3z) * u"A"
    end
    solution = 14.0u"A*m"

    # Package and run tests
    testable = TestableGeometry(integrand, ring, solution)
    runtests(testable)
end

@testitem "Meshes.Rope" setup=[Combinations] begin
    # Geometry
    a = Point(0, 0, 0)
    b = Point(1, 0, 0)
    c = Point(1, 1, 0)
    d = Point(1, 1, 1)
    rope = Rope(a, b, c, d)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        x, y, z = ustrip.((p.coords.x, p.coords.y, p.coords.z))
        (x + 2y + 3z) * u"A"
    end
    solution = 7.0u"A*m"

    # Package and run tests
    testable = TestableGeometry(integrand, rope, solution)
    runtests(testable)
end

@testitem "Meshes.Segment" setup=[Combinations] begin
    # Connect a line segment from the origin to an arbitrary point on the unit sphere
    φ, θ = (7pi / 6, pi / 3)  # Arbitrary spherical angles
    pt_a = Point(0, 0, 0)
    pt_b = Point(sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ))
    segment = Segment(pt_a, pt_b)

    # Integrand & Solution
    a, b = (7.1, 4.6)  # arbitrary constants > 0
    function integrand(p::Meshes.Point; a = a, b = b)
        r = ustrip(u"m", norm(to(p)))
        exp(r * log(a) + (1 - r) * log(b)) * u"A"
    end
    solution = ((a - b) / (log(a) - log(b))) * u"A*m"

    # Package and run tests
    testable = TestableGeometry(integrand, segment, solution)
    runtests(testable)
end

@testitem "Meshes.SimpleMesh" setup=[Combinations] begin
    # Geometry
    a = π
    points = [(0, 0), (a, 0), (0, a), (a, a), (0.25a, 0.5a), (0.75a, 0.5a)]
    tris = connect.([(1, 5, 3), (4, 6, 2)], Triangle)
    quads = connect.([(1, 2, 6, 5), (4, 3, 5, 6)], Quadrangle)
    mesh = SimpleMesh(points, [tris; quads])

    # Integrand & Solution
    function integrand(p::Meshes.Point; a = a)
        x₁, x₂ = ustrip.((to(p)))
        (√(a^2 - x₁^2) + √(a^2 - x₂^2)) * u"A"
    end
    solution = 2a * (π * a^2 / 4) * u"A*m^2"

    # Package and run tests
    testable = TestableGeometry(integrand, mesh, solution)
    runtests(testable; rtol = 1e-6)
end

@testitem "Meshes.Sphere 2D" setup=[Combinations] begin
    # Geometry
    origin = Point(0, 0)
    radius = 4.4
    sphere = Sphere(origin, radius)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        r = ustrip(u"m", norm(to(p)))
        exp(-r^2) * u"A"
    end
    solution = 2π * radius * exp(-radius^2) * u"A*m"

    # Package and run tests
    testable = TestableGeometry(integrand, sphere, solution)
    runtests(testable)
end

@testitem "Meshes.Sphere 3D" setup=[Combinations] begin
    # Geometry
    center = Point(1, 2, 3)
    radius = 4.4u"m"
    sphere = Sphere(center, radius)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        rθφ = convert(Spherical, Cartesian((p - center)...))
        r = ustrip(rθφ.r)
        θ = ustrip(rθφ.θ)
        φ = ustrip(rθφ.ϕ)
        (sin(φ)^2 + cos(θ)^2) * u"A"
    end
    solution = ((2π * radius^2) + (4π / 3 * radius^2)) * u"A"

    # Package and run tests
    testable = TestableGeometry(integrand, sphere, solution)
    runtests(testable)
end

@testitem "Meshes.Tetrahedron" setup=[Combinations] begin
    # Geometry
    a = Point(0, 0, 0)
    b = Point(1, 0, 0)
    c = Point(0, 1, 0)
    d = Point(0, 0, 1)
    tetrahedron = Tetrahedron(a, b, c, d)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        x, y, z = ustrip.(u"m", to(p))
        (x + 2y + 3z) * u"A"
    end
    solution = (1 // 4) * u"A*m^3"

    # Package and run tests
    testable = TestableGeometry(integrand, tetrahedron, solution)
    runtests(testable)
end

@testitem "Meshes.Torus" setup=[Combinations] begin
    # Geometry
    center = Point(0, 0, 0)
    ẑ = Vec(0, 0, 1)
    R = 3.5 # radius from axis-of-revolution to center of circle being revolved
    r = 1.2 # radius of circle being revolved
    torus = Torus(center, ẑ, R, r)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        x, y, z = ustrip.(u"m", to(p))
        (x^2 + y^2) * u"A"
    end
    solution = (2π^2 * r * R * (2R^2 + 3r^2)) * u"A*m^2"

    # Package and run tests
    testable = TestableGeometry(integrand, torus, solution)
    runtests(testable)
end

@testitem "Meshes.Triangle" setup=[Combinations] begin
    # Geometry
    a = Point(0, 0, 0)
    b = Point(1, 0, 0)
    c = Point(0, 1, 0)
    triangle = Triangle(a, b, c)

    # Integrand & Solution
    function integrand(p::Meshes.Point)
        x, y, z = ustrip.(u"m", to(p))
        (x + 2y + 3z) * u"A"
    end
    solution = (1 // 2) * u"A*m^2"

    # Package and run tests
    testable = TestableGeometry(integrand, triangle, solution)
    runtests(testable)
end

@testitem "Meshes.Wedge" setup=[Combinations] begin
    if pkgversion(Meshes) >= v"0.52.12"
        # Geometry
        a₀ = Point(0, 0, 0)
        b₀ = Point(1, 0, 0)
        c₀ = Point(0, 1, 0)
        a₁ = Point(0, 0, 1)
        b₁ = Point(1, 0, 1)
        c₁ = Point(0, 1, 1)
        wedge = Wedge(a₀, b₀, c₀, a₁, b₁, c₁)

        # Integrand & Solution
        function integrand(p::Meshes.Point)
            x, y, z = ustrip.(u"m", to(p))
            (x + 2y + 3z) * u"A"
        end
        solution = (5 // 4) * u"A*m^3"

        # Package and run tests
        testable = TestableGeometry(integrand, wedge, solution)
        runtests(testable)
    end
end
