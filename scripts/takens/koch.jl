# %%-----------Koch Curve---------------------

function pointskoch(points, maxk, α = sqrt(3)/2)
  Q = [0 -1; 1 0]
  for k = 1:maxk
    n = length(points)
    new_points = Vector{Float64}[]
    for i = 1:n-1
      p1, p2 = points[i], points[i+1]
      v = (p2 - p1) / 3
      q1 = p1 + v
      q2 = p1 + 1.5v + α * Q * v
      q3 = q1 + v
      append!(new_points, [p1, q1, q2, q3])
    end
    push!(new_points, points[end])
    points = new_points
  end
  return points
end
