################################################################################
#                            Common Methods
################################################################################

function lineintegral(
    f::F,
    ring::Meshes.Ring,
    settings::I
) where {F<:Function, I<:IntegrationAlgorithm}
    # Convert the Ring into Segments, sum the integrals of those 
    return sum(segment -> lineintegral(f, segment, settings), segments(ring))
end

function lineintegral(
    f::F,
    rope::Meshes.Rope,
    settings::I
) where {F<:Function, I<:IntegrationAlgorithm}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> lineintegral(f, segment, settings), segments(rope))
end


################################################################################
#                            Gauss-Legendre
################################################################################

"""
    lineintegral(f, curve::Meshes.BezierCurve, ::GaussLegendre; alg=Meshes.Horner())

Like [`lineintegral`](@ref) but integrates along the domain defined a `curve`.
By default this uses Horner's method to improve performance when parameterizing
the `curve` at the expense of a small loss of precision. Additional accuracy
can be obtained by specifying the use of DeCasteljau's algorithm instead with
`alg=Meshes.DeCasteljau()` but can come at a steep cost in memory allocations,
especially for curves with a large number of control points.
"""
function lineintegral(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V},
    settings::GaussLegendre;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function, Dim, T, V}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = gausslegendre(settings.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = 0.5x + 0.5
    point(x) = curve(t(x), alg)

    # Integrate f along the line and apply a domain-correction factor for [-1,1] ↦ [0, length]
    return 0.5 * length(curve) * sum(w .* f(point(x)) for (w,x) in zip(ws, xs))
end

function lineintegral(
    f::F,
    line::Meshes.Box{1,T},
    settings::GaussLegendre
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Box{1,T} is definitionally embedded in 1D-space
    _validate_integrand(f,1,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = gausslegendre(settings.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = 0.5x + 0.5
    point(x) = line(t(x))

    # Integrate f along the line and apply a domain-correction factor for [-1,1] ↦ [0, length]
    return 0.5 * length(line) * sum(w .* f(point(x)) for (w,x) in zip(ws, xs))
end

function lineintegral(
    f::F,
    circle::Meshes.Circle{T},
    settings::GaussLegendre
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Circle is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = gausslegendre(settings.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = 0.5x + 0.5
    point(x) = circle(t(x))

    # Integrate f along the circle's rim and apply a domain-correction
    #   factor for [-1,1] ↦ [0, circumference]
    return 0.5 * length(circle) * sum(w .* f(point(x)) for (w,x) in zip(ws, xs))
end

function lineintegral(
    f::F,
    line::Meshes.Line{Dim,T},
    settings::GaussLegendre
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = gausslegendre(settings.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = 0.5x + 0.5
    point(x) = line(t(x))

    # Integrate f along the line and apply a domain-correction factor for [-1,1] ↦ [0, length]
    return 0.5 * length(line) * sum(w .* f(point(x)) for (w,x) in zip(ws, xs))
end

function lineintegral(
    f::F,
    segment::Meshes.Segment{Dim,T},
    settings::GaussLegendre
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = gausslegendre(settings.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = 0.5x + 0.5
    point(x) = segment(t(x))

    # Integrate f along the line and apply a domain-correction factor for [-1,1] ↦ [0, length]
    return 0.5 * length(segment) * sum(w .* f(point(x)) for (w,x) in zip(ws, xs))
end

function lineintegral(
    f::F,
    circle::Meshes.Sphere{2,T},
    settings::GaussLegendre
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Sphere{2,T} is simply a circle in 2D-space
    _validate_integrand(f,2,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = gausslegendre(settings.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = 0.5x + 0.5
    point(x) = circle(t(x))

    # Integrate f along the circle's rim and apply a domain-correction
    #   factor for [-1,1] ↦ [0, circumference]
    return 0.5 * length(circle) * sum(w .* f(point(x)) for (w,x) in zip(ws, xs))
end


################################################################################
#                                Gauss-Kronrod
################################################################################

"""
    lineintegral(f, curve::BezierCurve, ::GaussKronrod; alg=Horner(), kws...)

Like [`lineintegral`](@ref) but integrates along the domain defined a `curve`.
By default this uses Horner's method to improve performance when parameterizing
the `curve` at the expense of a small loss of precision. Additional accuracy
can be obtained by specifying the use of DeCasteljau's algorithm instead with
`alg=Meshes.DeCasteljau()` but can come at a steep cost in memory allocations,
especially for curves with a large number of control points.
"""
function lineintegral(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V},
    settings::GaussKronrod;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function,Dim, T, V}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(curve)
    point(t) = curve(t, alg)
    return QuadGK.quadgk(t -> len * f(point(t)), 0, 1; settings.kwargs...)[1]
end

function lineintegral(
    f::F,
    line::Meshes.Box{1,T},
    settings::GaussKronrod
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Box is definitionally embedded in 1D-space
    _validate_integrand(f,1,T)

    len = length(line)
    point(t) = line(t)
    return QuadGK.quadgk(t -> len * f(point(t)), 0, 1; settings.kwargs...)[1]
end

function lineintegral(
    f::F,
    circle::Meshes.Circle{T},
    settings::GaussKronrod
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Circle is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    len = length(circle)
    point(t) = circle(t)
    return QuadGK.quadgk(t -> len * f(point(t)), 0, 1; settings.kwargs...)[1]
end

function lineintegral(
    f::F,
    line::Meshes.Line{Dim,T},
    settings::GaussKronrod
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(line)
    point(t) = line(t)
    return QuadGK.quadgk(t -> len * f(point(t)), 0, 1; settings.kwargs...)[1]
end

function lineintegral(
    f::F,
    segment::Meshes.Segment{Dim,T},
    settings::GaussKronrod
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(segment)
    point(t) = segment(t)
    return QuadGK.quadgk(t -> len * f(point(t)), 0, 1; settings.kwargs...)[1]
end

function lineintegral(
    f::F,
    circle::Meshes.Sphere{2,T},
    settings::GaussKronrod
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Sphere{2,T} is simply a circle in 2D-space
    _validate_integrand(f,2,T)

    len = length(circle)
    point(t) = circle(t)
    return QuadGK.quadgk(t -> len * f(point(t)), 0, 1; settings.kwargs...)[1]
end