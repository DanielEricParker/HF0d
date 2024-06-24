
function run_supermoire_model(;
        Umm=30,
        Umsm=30,
        Usmsm=0,
        W1=20,
        W2=1,
        δ=0.075,
        N=1000,
        repeats=12,
        μs = -150:0.25:150,
        save_all_species=false,
        verbose=false
    )
    Nf1,Nf2,Nf = 8,8,16
    hfm = HFmodel_supermoire(N; Umm=Umm, Umsm = Umsm, Usmsm = Usmsm, W1=W1, W2=W2, δ=δ);

    println("Running HF (0D). ")
    νsopt, Φsopt = run_HF(μs,hfm,repeats; verbose=verbose);

    #sort into mu order
    νsopt2 = stack(νsopt)
    νsopt2[1:Nf1,:] .= sort(νsopt2[1:Nf1,:],dims=1)
    νsopt2[Nf1+1:Nf,:] .= sort(νsopt2[Nf1+1:Nf,:],dims=1)

    νstot, κs = compressibility(μs, νsopt2; smooth=2, cutoff=1e-3)
    νstot .-= 4;

    dat = Dict() #h5open("data4/dat_W2_$(W2).h5","w") 
        dat["W1"] = W1
        dat["W2"] = W2
        dat["Umm"] = Umm
        dat["Umsm"] = Umsm
        dat["Usmsm"] = Usmsm
        dat["mus"] = collect(μs)
        if save_all_species
            dat["nus_all"] = νsopt2
        end
        dat["nus"] = νstot
        dat["kappas"] = κs
        dat["repeats"] = repeats
        dat["N"] = N
    return dat
end

function save_data_dict_to_dataframe!(datadict,filename)
    file = h5open(filename,"w")
    for (k,v) in datadict
        file[k] = v
    end
    close(file)
    return nothing
end

#call on collections of saved h5 files to process them to a dataframe
#there should be a more elegant way to do this...
function process_raw_data_files_to_dataframe(datafiles; verbose=true)
    dat = []
    for fn in datafiles
        @assert endswith(fn,".h5")
        file = h5open(fn, "r")
        verbose && println("Loading $(fn)...")
        datadict = (; (Symbol(k) => read(file,k) for k in keys(file))...)
        close(file)
        push!(dat,datadict)
    end
    df = DataFrame(dat);
    df[!,:W2] = [typeof(x)==String ? parse(Float64,x) : x for x in df[!,:W2]]
    df[!,:kappas2] = [max.(0.0,kappa) for kappa in df[!,:kappas]]
    return df
end


function generate_Wsm_scan_data(folder,W2s,μs = 0:0.25:150,verbose=true)
    println(Threads.nthreads())
    #embarrassingly parallel for loop we can multithread
    Threads.@threads for W2 in W2s
        try
            dat = run_supermoire_model(;
                    Umm=30,Umsm=30,Usmsm=0,
                    W1=20,W2=W2,δ=0.075,
                    N=1000, repeats=12,
                    μs = μs,
                    save_all_species=true,
                    verbose=verbose
                );
            filename=folder*"/dat_W2_$(W2).h5"
            save_data_dict_to_dataframe!(dat,filename)
        catch
            println("Failure at W2=$(W2)")
        end
    end
end

function generate_Um_scan_data(folder, Umms, μs = 0:0.25:150, verbose=true)
    println(Threads.nthreads())
    #embarrassingly parallel for loop we can multithread
    Threads.@threads for Um in Umms
        try
            dat = run_supermoire_model(;
                    Umm=Um,Umsm=Um,Usmsm=0,
                    W1=20,W2=5,δ=0.075,
                    N=1000, repeats=12,
                    μs = μs,
                    save_all_species=true,
                    verbose=verbose
                );
            filename=folder*"/dat_Um_$(Um).h5"
            save_data_dict_to_dataframe!(dat,filename)
        catch
            println("Failure at Um=$(Um)")
        end
    end
end