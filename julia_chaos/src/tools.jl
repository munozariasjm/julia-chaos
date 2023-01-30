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
function get_line_params(x1, y1, x2, y2)
    slope = (y1 - y2) / (x1 - x2)
    intercept = y1 - slope * x1
    return slope, intercept
end
function euclidean_distance(point1, point2)
    return √((point1[1] - point2[1])^2 + (point1[2] - point2[2])^2)
end

function solve_quadratic_equation(A, B, C)
    x1 = (-B + √(B^2 - 4 * A * C)) / (2 * A)
    x2 = (-B - √(B^2 - 4 * A * C)) / (2 * A)
    return [x1, x2]
end
function get_datapoint(pivot, measure, length, direction="inner")
    if pivot[1] == measure[1]
        y1 = pivot[2] + length
        y2 = pivot[2] - length
        x1 = pivot[1]
        x2 = pivot[1]
    else
        slope, intercept = get_line_params(pivot[1], pivot[2],
            measure[1], measure[2],)
        A = 1
        B = -2 * pivot[1]
        C = pivot[1]^2 - length^2 / (slope^2 + 1)
        x1, x2 = solve_quadratic_equation(A, B, C)
        y1 = slope * x1 + intercept
        y2 = slope * x2 + intercept
    end
    if direction == "inner"
        if euclidean_distance((x1, y1), measure) < euclidean_distance((x2, y2), measure)
            datapoint = (x1, y1)
        else
            datapoint = (x2, y2)
        end
    else
        if euclidean_distance((x1, y1), measure) > euclidean_distance((x2, y2), measure)
            datapoint = (x1, y1)
        else
            datapoint = (x2, y2)
        end
    end
    return datapoint
end

export analytic_solve, is_stalbe, meshgrid, rungekutta4, get_datapoint, euclidean_distance, solve_quadratic_equation

end