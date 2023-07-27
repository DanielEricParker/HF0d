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

function HFmodel_supermoire(N; U0 = 30,U1=15, δ=0.075)
    Nf1,Nf2 = 8,8
    Nf = Nf1+Nf2

    ϵs = Array{Array{Float64,1},1}(undef,Nf)
    ρs = Array{Array{Float64,1},1}(undef,Nf)
    νs = Array{Array{Float64,1},1}(undef,Nf)
    Ekins = Array{Array{Float64,1},1}(undef,Nf)

    n1total = (1-δ*Nf2/Nf1) #2 bands
    n2total = δ     #2 band
    totaldensities = vcat(fill(n1total,Nf1),fill(n2total,Nf2))

    W1,W2 = 20, 5
    ϵs1 = 0:(W1/N):W1
    ϵs2 = -5:(W1/N):5
    ρvanHove1(ϵ) =  min(1,ρvanHove(ϵ,W1,W1,0.7*W1))
    ρvanHove2(ϵ) =  min(2,ρvanHove(ϵ,W2,W2,0.7*W2))

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
    end
    for n in 1:div(Nf1,2)
        Ekins[n] .-= Ekins[n][end]
    end
    for α in Nf1+1:Nf
        Ekins[α] = cumualtivetrapezoidintegration(ϵs[α], abs.(ϵs[α]) .* ρs[α])
        Ekins[α] .-= Ekins[α][div(length(Ekins[α]),2)+1];
    end
    Ekin_ns = [linear_interpolation(νs[α],Ekins[α], extrapolation_bc=Line()) for α in 1:Nf];

    Umat = kron([U0 U0; U0 U1], fill(1.0,Nf1,Nf1)-I)

    return HFModel(Nf, Umat, ϵs, ρs, νs, Ekins, Ekin_ns);
end


Ekinetic(νs,Es) = sum(E(ν) for (ν, E) in zip(νs,Es))
function grandpotential(νs, μ, hfm)
    @assert length(νs) == hfm.Nf
    E = Ekinetic(νs,hfm.Ekin_ns) #kinetic 
    E += dot(νs,hfm.Umat,νs)/2    #potential
    E -= μ * sum(νs)    #chemical potential
    return E
end

function run_HF(μs, hfm, repeats = 10, verbose=true)
    νsopt = Array{Array{Float64,1},1}([])
    Φsopt = Array{Float64,1}([])

    lowerbounds = [νs[1] for νs in hfm.νs]
    upperbounds = [νs[end] for νs in hfm.νs]
    random_ναs() = [rand(Uniform(lowerbounds[α],upperbounds[α])) for α in 1:hfm.Nf]
    for (n,μ) in enumerate(μs)
        verbose && n % 10 == 0 && println("μ = $μ ($n/$(length(μs)))")
        F(μs) = grandpotential(μs,μ,hfm)
        #repeat several times to find global minimum
        results = [optimize(F,lowerbounds,upperbounds,random_ναs()) for _ in 1:repeats]
        opt = argmin(res -> Optim.minimum(res), results) #select global minimum 
        push!(νsopt, Optim.minimizer(opt))
        push!(Φsopt, Optim.minimum(opt))
    end
    return νsopt, Φsopt
end

