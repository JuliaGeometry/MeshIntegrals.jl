using BenchmarkTools
using LinearAlgebra
using Meshes
using MeshIntegrals
using StaticArrays
using Unitful
import Enzyme

const SUITE = BenchmarkGroup()

############################################################################################
#                                      Integrals
############################################################################################

# Like fill but returns an SVector to eliminate integrand function allocations
staticfill(val, N) = SVector(ntuple(Returns(val), N)...)

integrands = (
    (name = "Scalar", f = p -> norm(to(p))),
    (name = "Vector", f = p -> staticfill(norm(to(p)), 3))
)
rules = (
    (name = "GaussLegendre", rule = GaussLegendre(100)),
    (name = "GaussKronrod", rule = GaussKronrod()),
    (name = "HAdaptiveCubature", rule = HAdaptiveCubature())
)
geometries = (
    (name = "Segment", item = Segment(Point(0, 0, 0), Point(1, 1, 1))), # 1D
    (name = "Sphere", item = Sphere(Point(0, 0), 1.0)), # 2D
    (name = "Sphere", item = Sphere(Point(0, 0, 0), 1.0)) # 3D
)

SUITE["Integrals"] = let s = BenchmarkGroup()
    for (int, rule, geometry) in Iterators.product(integrands, rules, geometries)
        # Skip unsupported applications of GaussKronrod
        if (Meshes.paramdim(geometry.item) > 1) && (rule.rule == GaussKronrod())
            continue
        end

        # Construct benchmark and add it to test suite
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
    g = (
        bezier = BezierCurve([Point(t, sin(t), 0) for t in range(-pi, pi, length = 361)]),
        line = Line(Point(0, 0, 0), Point(1, 1, 1)),
        plane = Plane(Point(0, 0, 0), Vec(0, 0, 1)),
        ray = Ray(Point(0, 0, 0), Vec(0, 0, 1)),
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
    rule_gl = GaussLegendre(100),
    rule_gk = GaussKronrod(),
    rule_hc = HAdaptiveCubature()
)

SUITE["Specializations/Scalar GaussLegendre"] = let s = BenchmarkGroup()
    s["BezierCurve"] = @benchmarkable integral($spec.f, $spec.g.bezier, $spec.rule_gl)
    s["Line"] = @benchmarkable integral($spec.f_exp, $spec.g.line, $spec.rule_gl)
    s["Plane"] = @benchmarkable integral($spec.f_exp, $spec.g.plane, $spec.rule_gl)
    s["Ray"] = @benchmarkable integral($spec.f_exp, $spec.g.ray, $spec.rule_gl)
    s["Rope"] = @benchmarkable integral($spec.f, $spec.g.rope, $spec.rule_gl)
    s["Triangle"] = @benchmarkable integral($spec.f, $spec.g.triangle, $spec.rule_gl)
    s["Tetrahedron"] = @benchmarkable integral($spec.f, $spec.g.tetrahedron, $spec.rule_gl)
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
