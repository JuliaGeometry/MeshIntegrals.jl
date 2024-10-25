using BenchmarkTools
using LinearAlgebra
using Meshes
using MeshIntegrals
using Unitful

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
    (name = "Meshes.Segment", item = Segment(Point(0, 0, 0), Point(1, 1, 1))),
    (name = "Meshes.Sphere", item = Sphere(Point(0, 0, 0), 1.0))
)

SUITE["Integrals"] = let s = BenchmarkGroup()
    for (int, rule, geometry) in Iterators.product(integrands, rules, geometries)
        n1, n2, N = geometry.name, "$(int.name) $(rule.name)", rule.N
        s[n1, n2] = @benchmarkable integral($int.f, $geometry.item, $rule.rule) evals=N
    end
    s
end

############################################################################################
#                                      Differentials
############################################################################################

sphere = geometries[2].item
differential = MeshIntegrals.differential

SUITE["Differentials"] = let s = BenchmarkGroup()
    s["Jacobian"] = @benchmarkable jacobian($sphere, $(0.5, 0.5)) evals=1000
    s["Differential"] = @benchmarkable differential($sphere, $(0.5, 0.5)) evals=1000
    s
end

#tune!(SUITE)
#run(SUITE, verbose=true)
