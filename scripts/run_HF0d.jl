using LinearAlgebra
using HDF5
using Interpolations

push!(LOAD_PATH,"../")
using HF0d


function runHF(;Umm=30,Umsm=30,Usmsm=0,W1=20,W2=20,δ=0.075,N=1000,repeats=12)
    Nf1,Nf2,Nf = 8,8,16
    hfm = HF0d.HFmodel_supermoire(N; Umm=Umm, Umsm = Umsm, Usmsm = Usmsm, W1=W1, W2=W2, δ=δ);
    # μs = range(120,250,length=NN)
    μs = 120:0.25:240 #120:2:250

    println("Running HF (0D). ")
    νsopt, Φsopt = run_HF(μs,hfm,repeats; verbose=false);

    νsopt2 = stack(νsopt)
    νsopt2[1:Nf1,:] .= sort(νsopt2[1:Nf1,:],dims=1)
    νsopt2[Nf1+1:Nf,:] .= sort(νsopt2[Nf1+1:Nf,:],dims=1)

    νstot, κs = compressibility(μs, νsopt2; smooth=2, cutoff=1e-3)
    νstot .-= 4;

    file = h5open("data4/dat_W2_$(W2).h5","w") 
        file["W1"] = W1
        file["W2"] = W2
        file["Umm"] = Umm
        file["Umsm"] = Umsm
        file["Usmsm"] = Usmsm
        file["mus"] = collect(μs)
        file["nus"] = νstot
        file["kappas"] = κs
        # file["nus_interp"] = collect(νs_interp)
        # file["kappas_interp"] = κs_interp
        file["repeats"] = repeats
        file["N"] = N
        file["NN"] = NN
    close(file)
    return nothing
end

# runHF()

function main()
    println(Threads.nthreads())
    W2s = 16:1:24
    Threads.@threads for W2 in W2s
        println(Threads.threadid())
        try
            runHF(W2=W2)
        catch
            println("Failure at W2=$(W2)")
        end
    end
end
main()