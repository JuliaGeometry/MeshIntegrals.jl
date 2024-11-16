using BenchmarkTools
using LinearAlgebra
using Meshes
using MeshIntegrals

const SUITE = BenchmarkGroup()

############################################################################################
#                                      Integrals
############################################################################################

integrands = (
    (name = "Scalar", f = p -> norm(to(p))),
    (name = "Vector", f = p -> fill(norm(to(p)), 3))
)
rules = (
    (name = "GaussLegendre", rule = GaussLegendre(100), N = 100),
    (name = "GaussKronrod", rule = GaussKronrod(), N = 100),
    (name = "HAdaptiveCubature", rule = HAdaptiveCubature(), N = 500)
)
geometries = (
    (name = "Segment", item = Segment(Point(0, 0, 0), Point(1, 1, 1))),
    (name = "Sphere", item = Sphere(Point(0, 0, 0), 1.0))
)

SUITE["Integrals"] = let s = BenchmarkGroup()
    for (int, rule, geometry) in Iterators.product(integrands, rules, geometries)
        n1, n2, N = geometry.name, "$(int.name) $(rule.name)", rule.N
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
    g = (
        bezier = BezierCurve([Point(t, sin(t), 0) for t in range(-pi, pi, length = 361)]),
        line = Line(Point(0, 0, 0), Point(1, 1, 1)),
        plane = Plane(Point(0, 0, 0), Vec(0, 0, 1)),
        ray = Ray(Point(0, 0, 0), Vec(0, 0, 1)),
        triangle = Triangle(Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1)),
        tetrahedron = let
            a = Point(0, 3, 0)
            b = Point(-7, 0, 0)
            c = Point(8, 0, 0)
            ẑ = Vec(0, 0, 1)
            Tetrahedron(a, b, c, a + ẑ)
        end
    ),
    rule = GaussLegendre(100)
)

SUITE["Specializations/Scalar GaussLegendre"] = let s = BenchmarkGroup()
    s["BezierCurve"] = @benchmarkable integral($spec.f, $spec.g.bezier, $spec.rule)
    s["Line"] = @benchmarkable integral($spec.f_exp, $spec.g.line, $spec.rule)
    s["Plane"] = @benchmarkable integral($spec.f_exp, $spec.g.plane, $spec.rule)
    s["Ray"] = @benchmarkable integral($spec.f_exp, $spec.g.ray, $spec.rule)
    s["Triangle"] = @benchmarkable integral($spec.f, $spec.g.triangle, $spec.rule)
    s["Tetrahedron"] = @benchmarkable integral($spec.f, $spec.g.tetrahedron, $spec.rule)
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

#tune!(SUITE)
#run(SUITE, verbose=true)
