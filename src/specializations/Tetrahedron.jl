################################################################################
#                  Specialized Methods for Tetrahedron
################################################################################

function integral(
    f::F,
    tetrahedron::Meshes.Tetrahedron,
    settings::GaussLegendre,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    error("Integrating a Tetrahedron with GaussLegendre not supported.")
end

function integral(
    f::F,
    tetrahedron::Meshes.Tetrahedron,
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    inner∫₂(v,w) = QuadGK.quadgk(u -> f(tetrahedron(u,v,w)), FP(0), FP(1-v-w); settings.kwargs...)[1]
    inner∫₁(w) = QuadGK.quadgk(v -> inner∫₂(v,w), FP(0), FP(1-w); settings.kwargs...)[1]
    outer∫ = QuadGK.quadgk(w -> inner∫₁(w), FP(0), FP(1); settings.kwargs...)[1]

    # Apply barycentric domain correction (volume: 1/6 → actual)
    return 6 * volume(tetrahedron) * outer∫
end

function integral(
    f::F,
    tetrahedron::Meshes.Tetrahedron,
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64,
) where {F<:Function, T<:AbstractFloat}
    error("Integrating a Tetrahedron with HAdaptiveCubature not supported.")
end
