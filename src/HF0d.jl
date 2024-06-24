module HF0d

using LinearAlgebra
using StaticArrays
using Interpolations
using Distributions
using Optim
using ImageFiltering
using DataFrames
using HDF5

include("utilities.jl")
include("HFmodel.jl")


export trapezoidintegration,cumualtivetrapezoidintegration
export HFModel, grandpotential, ∇grandpotential!, run_HF,run_HF2
export HFModel_linearDOS, HFModel_vanHoveDOS, HFmodel_supermoire
export compressibility

include("run_and_process_data.jl")
export run_supermoire_model, save_data_dict_to_dataframe!, process_raw_data_files_to_dataframe
export generate_Wsm_scan_data, generate_Um_scan_data

end # module trilayer_SET_model
