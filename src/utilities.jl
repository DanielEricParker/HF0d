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


# function _rhobar(νs,U,ρ_νs, lbs, ubs)
#     work = zero(Float64)
#     for (ν,ρ_ν,lb,ub) in zip(νs,ρ_νs,lbs,ubs)
#         ρα = ifelse(lb*0.99999 < ν < ub * 0.99999, ρ_ν(ν), 0) #dμ/dν vanishes in this case
#         work += ρα/(1-ρα*U)
#     end
#     return work
# end

# # Zondinger et al method
# function compressibility(νsopt,hfm)
#     #assumes all values of U are the same
#     U = hfm.Umat[1,2]
#     lb = [νs[1] for νs in hfm.νs]
#     ub = [νs[end] for νs in hfm.νs]
#     ρs_ν = [linear_interpolation(ν,ρ, extrapolation_bc=Line()) for (ν,ρ) in zip(hfm.νs,hfm.ρs)]
#     ρbars = [_rhobar(νs_μ,U,ρs_ν,lb,ub) for νs_μ in νsopt]
#     κs = @. ρbars /(1+ρbars*U)
#     return κs
# end

#analytically compute dn/dμ
#doesn't work for the two-species case for unclear reasons

# function compressibility2(νsopt,hfm)
#     κs = zeros(Float64, length(νsopt))
#     Ut = copy(hfm.Umat)
#     v = ones(Float64,hfm.Nf)
#     @assert norm(diag(Ut)) ≈ zero(Float64)
#     lbs = [νs[1] for νs in hfm.νs]
#     ubs = [νs[end] for νs in hfm.νs]
#     ρs_ν = [linear_interpolation(ν,ρ, extrapolation_bc=Line()) for (ν,ρ) in zip(hfm.νs,hfm.ρs)]
#     for n in eachindex(κs)
#         Ut[:] = hfm.Umat
#         κs[n] = _comp_helper(Ut,v, νsopt[n], ρs_ν,lbs,ubs)
#     end
#     return κs
# end

# function _comp_helper(Ut,v,νs,ρs,lbs,ubs)
#     mask = zeros(Bool,hfm.Nf)
#     for α in axes(Ut,1)
#         if lbs[α]*0.9999 < νs[α] < ubs[α] * 0.9999
#             Ut[α,α] = 1/ρs[α](νs[α])
#             mask[α] = true #dμ/dν non-vanishing
#         end
#     end
#     v2 = Ut[mask,mask] \ v[mask]
#     return sum(v2)
# end

#just take the discrete differences with a bit of smoothing --- much easier
function compressibility(μs, νsopt; smooth = 2, cutoff=1e-3)
    ker = ImageFiltering.Kernel.gaussian((smooth,))
    νsopt = mapslices(n -> imfilter(n,ker), νsopt, dims=2);
    νstot = dropdims(sum(νsopt,dims=1),dims=1)
    κs = (diff(νstot) ./ diff(μs)) .+ cutoff
    κs .*= (νstot[end]-νstot[1])/trapezoidintegration(μs[1:end-1], κs) 
    return νstot[1:end-1], κs
end