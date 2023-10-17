"""
Represents a 0-dimensional Hartree Fock Model
"""
struct HFModel{T}
    Nf      :: Int64      #number of flavors
    Umat    :: Array{Float64,2}  #interaction array
    ϵs      :: Array{Array{Float64,1},1} 
    ρs      :: Array{Array{Float64,1},1}
    νs      :: Array{Array{Float64,1},1}
    Ekins   :: Array{Array{Float64,1},1}
    Ekin_ns :: T        #function or interpolation E_{kinetic}(μ_α)
    ν0      :: Array{Float64,1}
    dν0     :: Array{Int64,1}
end

raw"""Given the number of flavors `Nf`, each with energies `ϵ` density of states `ρ` so that
``\int \rho(\epsilon) d\epsilon = n_{tot}``, the total density of that flavor, and the
interaction matrix ``U_{\alpha\beta}``, creates the corresponding HarteeFockModel
"""
function HFModel(Nf,ϵs,ρs,totaldensities, Umat)  
    @assert length(ϵs) == Nf
    νs = Array{Array{Float64,1},1}(undef,Nf)
    Ekins = Array{Array{Float64,1},1}(undef,Nf)
    for α in 1:Nf
        @assert trapezoidintegration(ϵs[α],ρs[α]) ≈ totaldensities[α]
        νs[α] = zeromiddle!(cumualtivetrapezoidintegration(ϵs[α], ρs[α]))
        Ekins[α] = zeromiddle!(cumualtivetrapezoidintegration(ϵs[α], ϵs[α] .* ρs[α]))
    end
    Ekin_ns = [linear_interpolation(νs[α],Ekins[α], extrapolation_bc=Line()) for α in 1:Nf]
    return HFModel(Nf, Umat, ϵs, ρs, νs, Ekins, Ekin_ns)
end

# various DOS models 
ρlinear(ϵ,W) = 2*abs(ϵ)/W^2
ρvanHove(ϵ, W, ϵ0, ϵv) =  (abs(ϵ)/W^2)*abs(ϵ0/(ϵv - abs(ϵ)))^(1/2)
ρvanHoveSM(ϵ, W, ϵ0, ϵv) =  (abs(ϵ)/W^2)*abs(ϵ0/(ϵv - abs(ϵ)))^(1)


"""
Constructs the Hartree-Fock Model with linear density of states from Supplement SI11 of
`Cascade of phase transitions and Dirac revivals in magic-angle graphene` by Zondinger et al.
"""
function HFModel_linearDOS(Nf,N,W,U)
    Umat = U * (fill(1.0,(Nf,Nf)) - I)

    ϵs = [collect(-W : (W/N) : W) for α in 1:Nf]
    ρs = [map(e -> ρlinear(e,W),ϵs[α]) for α in 1:Nf]
    total_densities = trapezoidintegration.(ϵs,ρs)

    return HFModel(Nf,ϵs,ρs,total_densities,Umat)
end

"""
Constructs the Hartree-Fock Model with a van Hove singularity from Supplement SI12 of
`Cascade of phase transitions and Dirac revivals in magic-angle graphene` by Zondinger et al.
"""
function HFModel_vanHoveDOS(Nf,N,W,U)
    Umat = U * (fill(1.0,(Nf,Nf)) - I)

    ϵs = [collect(-W : (W/N) : W) for α in 1:Nf]
    ρs = map(ϵs) do ϵα
        ρα = min.(map(e->ρvanHove(e,W,W,0.7*W),ϵα), 5)
        ρα .*= 2/trapezoidintegration(ϵα,ρα)
    end
    total_densities = trapezoidintegration.(ϵs,ρs)

    return HFModel(Nf,ϵs,ρs,total_densities,Umat)
end


function HFmodel_supermoire(N; Umm = 30, Usmsm=15, Umsm = 30, W1=20, W2 = 20, δ=0.075)
    Nf1,Nf2 = 8,8
    Nf = Nf1+Nf2
    ra1m = 1:div(Nf1,2)
    ra1p = div(Nf1,2)+1:Nf1
    ra2m = Nf1+1:Nf1+div(Nf2,2)
    ra2p = Nf1+div(Nf2,2)+1:Nf1+Nf2

    ϵs = Array{Array{Float64,1},1}(undef,Nf)
    ρs = Array{Array{Float64,1},1}(undef,Nf)
    νs = Array{Array{Float64,1},1}(undef,Nf)
    Ekins = Array{Array{Float64,1},1}(undef,Nf)

    n1total = (1-δ*Nf2/Nf1) #2 bands
    n2total = δ     #2 band
    totaldensities = vcat(fill(n1total,Nf1),fill(n2total,Nf2))
    ν0 = deepcopy(totaldensities)
    dν0 = zeros(Int64,length(ν0))
    # W1,W2 = 20, 15
    ϵs1 = 0:(W1/N):W1
    ϵs2 = 0:(W2/N):W2
    ρvanHove1(ϵ) =  min(1,ρvanHove(ϵ,W1,W1,0.7*W1))
    ρvanHove2(ϵ) =  min(10,ρvanHoveSM(ϵ,W2,W2,0.7*W2))

    #initialize bands 
    for (n1,n2) in zip(4:-1:1,5:8)
        ρs1 = ρvanHove1.(ϵs1)
        ρs1 .*= n1total/trapezoidintegration(ϵs1,ρs1)
        #hole bands 
        ϵs[n1] = collect(-reverse(ϵs1))
        ρs[n1] = reverse(ρs1)
        #electron bands
        ϵs[n2] = ϵs1
        ρs[n2] = ρs1

        ρs2 = ρvanHove2.(ϵs2)
        ρs2 .*= n2total/trapezoidintegration(ϵs2,ρs2)
        #hole bands
        ϵs[n1+Nf1] = collect(-reverse(ϵs2))
        ρs[n1+Nf1] = reverse(ρs2)
        #electron bands
        ϵs[n2+Nf1] = ϵs2
        ρs[n2+Nf1] = ρs2
    end
    
    for α in 1:Nf
        @assert trapezoidintegration(ϵs[α],ρs[α]) ≈ totaldensities[α]
        νs[α] = cumualtivetrapezoidintegration(ϵs[α], ρs[α]) #uncentered
        Ekins[α] = cumualtivetrapezoidintegration(ϵs[α], (ϵs[α]) .* ρs[α]) #uncentered
        if α in ra1m || α in ra2m
            Ekins[α] .-= Ekins[α][end]
            dν0[α] = -1
        else
            ν0[α] =  0
            dν0[α] = 1
        end
    end


    Ekin_ns = [linear_interpolation(νs[α],Ekins[α], extrapolation_bc=Line()) for α in 1:Nf];

    Umat = kron([Umm Umsm; Umsm Usmsm], fill(1.0,Nf1,Nf1)-I)

    return HFModel(Nf, Umat, ϵs, ρs, νs, Ekins, Ekin_ns,ν0,dν0);
end


Ekinetic(νs,Es) = sum(E(ν) for (ν, E) in zip(νs,Es))
function grandpotential(νs, μ, hfm)
    @assert length(νs) == hfm.Nf
    E = Ekinetic(νs,hfm.Ekin_ns) #kinetic
    ns = hfm.dν0 .* (νs - hfm.ν0)
    E += dot(ns,hfm.Umat,ns)/2    #potential
    E -= μ * sum(νs)    #chemical potential
    return E
end

function ∇grandpotential!(∇Φ,νs,μ,hfm)
    for i in eachindex(∇Φ)
        Interpolations.gradient(hfm.Ekin_ns[i],νs[i])[1]
        ∇Φ[i] = Interpolations.gradient(hfm.Ekin_ns[i],νs[i])[1]
        ∇Φ[i] -= μ
    end
    ns = hfm.dν0 .* (νs - hfm.ν0)
    ∇Φ .+= hfm.dν0 .* (hfm.Umat * ns)
    return ∇Φ
end

function run_HF(μs, hfm, repeats = 12; verbose=true)
    νsopt = Array{Array{Float64,1},1}([])
    Φsopt = Array{Float64,1}([])

    lowerbounds = [νs[1] for νs in hfm.νs]
    upperbounds = [νs[end] for νs in hfm.νs]
    random_ναs() = [rand(Uniform(lowerbounds[α],upperbounds[α])) for α in 1:hfm.Nf]
    for (n,μ) in enumerate(μs)
        verbose && n % 10 == 0 && println("μ = $μ ($n/$(length(μs)))")
        F(μs) = grandpotential(μs,μ,hfm)
        dF!(buffer, μs) = ∇grandpotential!(buffer,μs,μ,hfm)
        
        results = []
        for m in 1:repeats
            if m==1 & n >1
                ν0 = νsopt[n-1]
            else
                ν0 = random_ναs()
            end
            opt = optimize(F, dF!, lowerbounds, upperbounds,ν0 , Fminbox(GradientDescent()))
            if !isnan(Optim.minimum(opt))
                push!(results,opt)
            end
        end
        opt = argmin(res -> Optim.minimum(res),results) #select global minimum 
        push!(νsopt, Optim.minimizer(opt))
        push!(Φsopt, Optim.minimum(opt))
    end
    return νsopt, Φsopt
end