function trapezoidintegration(X,Y)
    @assert length(X) == length(Y)
    S = promote_type(eltype(X),eltype(Y))
    out = zero(S)
    for i in eachindex(X,Y)
        i == 1 && continue
        out += 0.5 * (X[i]-X[i-1]) *(Y[i]+Y[i-1])
    end
    return out
end

function cumualtivetrapezoidintegration(X,Y)
    @assert length(X) == length(Y)
    S = promote_type(eltype(X),eltype(Y))
    out = zeros(S,length(X))
    for i in eachindex(X,Y)
        i == 1 && continue
        out[i] = out[i-1] + 0.5* (X[i]-X[i-1]) *(Y[i]+Y[i-1])
    end
    return out
end

zeromiddle!(arr) = (arr .-= arr[div(length(arr),2)+1])

#take the discrete differences with a bit of smoothing due to finite resolution
function compressibility(μs, νsopt; smooth = 2, cutoff=1e-3)
    ker = ImageFiltering.Kernel.gaussian((smooth,))
    νsopt = mapslices(n -> imfilter(n,ker), νsopt, dims=2);
    νstot = dropdims(sum(νsopt,dims=1),dims=1)
    κs = (diff(νstot) ./ diff(μs)) .+ cutoff
    κs .*= (νstot[end]-νstot[1])/trapezoidintegration(μs[1:end-1], κs) 
    return νstot[1:end-1], κs
end