using BenchmarkTools
using LinearAlgebra
using Meshes
using MeshIntegrals
using Unitful
import Enzyme

const SUITE = BenchmarkGroup()

############################################################################################
#                                      Integrals
############################################################################################

integrands = (
    (name = "Scalar", f = p -> norm(to(p))),
    (name = "Vector", f = p -> fill(norm(to(p)), 3))
)
rules = (
    (name = "GaussLegendre", rule = GaussLegendre(100)),
    (name = "GaussKronrod", rule = GaussKronrod()),
    (name = "HAdaptiveCubature", rule = HAdaptiveCubature())
)
geometries = (
    (name = "Segment", item = Segment(Point(0, 0, 0), Point(1, 1, 1))),
    (name = "Sphere", item = Sphere(Point(0, 0, 0), 1.0))
)

SUITE["Integrals"] = let s = BenchmarkGroup()
    for (int, rule, geometry) in Iterators.product(integrands, rules, geometries)
        n1 = geometry.name
        n2 = "$(int.name) $(rule.name)"
        s[n1][n2] = @benchmarkable integral($int.f, $geometry.item, $rule.rule)
    end
    s
end

############################################################################################
#                                    Specializations
############################################################################################

spec = (
    f = p -> norm(to(p)),
    f_exp = p::Point -> exp(-norm(to(p))^2 / u"m^2"),
    geometries = (
        bezier = BezierCurve([Point(t, sin(t), 0) for t in range(-pi, pi, length = 361)]),
        rope = Rope([Point(t, t, t) for t in 1:32]...),
        triangle = Triangle(Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1)),
        tetrahedron = let
            a = Point(0, 3, 0)
            b = Point(-7, 0, 0)
            c = Point(8, 0, 0)
            ẑ = Vec(0, 0, 1)
            Tetrahedron(a, b, c, a + ẑ)
        end
    ),
    geometries_exp = (
        line = Line(Point(0, 0, 0), Point(1, 1, 1)),
        plane = Plane(Point(0, 0, 0), Vec(0, 0, 1)),
        ray = Ray(Point(0, 0, 0), Vec(0, 0, 1))
    ),
    rules = (
            ( name = "GaussLegendre", rule = GaussLegendre(100) ),
            ( name = "HAdaptiveCubature", rule = HAdaptiveCubature() )
    )
)

SUITE["Specializations"] = let s = BenchmarkGroup()
    s[]
    for r in spec.rules, geometry in spec.geometries
        geometry_name = nameof(typeof(geometry))
        s[geometry_name, r.name] = @benchmarkable integral($spec.f, geometry, r.rule)
    end
    for r in spec.rules, geometry in spec.geometries_exp
        geometry_name = nameof(typeof(geometry))
        s[geometry_name, r.name] = @benchmarkable integral($spec.f_exp, geometry, r.rule)
    end
    s
end

############################################################################################
#                                      Differentials
############################################################################################

sphere = Sphere(Point(0, 0, 0), 1.0)
differential = MeshIntegrals.differential

SUITE["Differentials"] = let s = BenchmarkGroup()
    s["Jacobian"] = @benchmarkable jacobian($sphere, $(0.5, 0.5)) evals=10
    s["Differential"] = @benchmarkable differential($sphere, $(0.5, 0.5)) evals=10
    s
end

############################################################################################
#                                   Integration Rules
###########################################################################################

SUITE["Rules"] = let s = BenchmarkGroup()
    s["GaussLegendre"] = @benchmarkable GaussLegendre($1000)
    s
end

#tune!(SUITE)
#run(SUITE, verbose=true)
