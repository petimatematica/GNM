function points_on_circle(C, r, np)
    alpha = range(0.0, stop=2π, length=np+1)
    points = Matrix{Float64}(undef, np, 2)
    for i in 1:np
        x = C[1] + r * cos(alpha[i])
        y = C[2] + r * sin(alpha[i])
        points[i, :] = [x, y]
    end

    return points
end

