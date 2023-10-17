module HF0d

using LinearAlgebra
using StaticArrays
using Interpolations
using Distributions
using Optim
using ImageFiltering

include("utilities.jl")
include("HFmodel.jl")

export trapezoidintegration,cumualtivetrapezoidintegration
export HFModel, grandpotential, ∇grandpotential!, run_HF,run_HF2
export HFModel_linearDOS, HFModel_vanHoveDOS, HFmodel_supermoire
export compressibility

end # module trilayer_SET_model
