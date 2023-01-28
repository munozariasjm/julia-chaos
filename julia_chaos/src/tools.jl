module ChaosTools
using Plots
"""
analytic_solve: Analitically solve any one dimensional flow system
"""
function is_stalbe(x::Float64, f::Function; δ=0.1)
    # Check if the left goes up and the right goes to the left
    δ = abs(δ * x + 1e-10)
    if f(x - δ) > f(x) && f(x + δ) < f(x)
        return true
    end
    return false
end

meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

function analytic_solve(f::Function, xs::AbstractArray; n::Int=1000, δ::Float64=0.00001, title="")
    # Plot the
    max_x = maximum(xs)
    min_x = minimum(xs)
    ys = f.(xs)
    plot(xs, ys, label="f(x)", title=title, xlim=(max_x, min_x))
    # Plot the fixed points
    positive_mask = ys .> 0
    # Plot the positive values
    scatter!(xs[positive_mask], ys[positive_mask] .* 0, label="positive", markershape=:rtriangle)
    # Keep only the negative values
    negative_mask = .!positive_mask
    # Plot the negative values
    scatter!(xs[negative_mask], ys[negative_mask] .* 0, label="negative", markershape=:ltriangle)
    # Scatter the points where the function is zero
    closer_to_zero = abs.(f.(range(minimum(xs), maximum(xs), n))) .< δ
    interscet = range(minimum(xs), maximum(xs), n)[closer_to_zero]
    scatter!(interscet, f.(interscet), label="Fixed points", markershape=:circle)
    # Check if the fixed points are stable
    stable = is_stalbe.(interscet, f)
    scatter!(interscet[stable], f.(interscet[stable]), label="Sink", markershape=:circle)
end
function rungekutta4(f::Function, y0::Vector, t::Vector)
    y = zeros((length(t), length(y0)))
    y[1, :] = y0
    for i in eachindex(t[1:end-1])
        h = t[i+1] - t[i]
        k1 = f(y[i, :], t[i])
        k2 = f(y[i, :] + k1 * h / 2, t[i] + h / 2)
        k3 = f(y[i, :] + k2 * h / 2, t[i] + h / 2)
        k4 = f(y[i, :] + k3 * h, t[i] + h)
        y[i+1, :] = y[i, :] + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    end
    return y
end

export analytic_solve, is_stalbe, meshgrid, rungekutta4

end