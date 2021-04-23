using Pkg
using DataFrames, DataStructures, CSV
# using Plots,StatsPlots
using Gadfly, Cairo, Fontconfig, ColorSchemes, Colors
using Random, StatsBase, Statistics
using Contour


# replace this with the directory for your own machine
home = "/Users/alexandermeyer/Dropbox/Research/Project02_VM/Code/upwelling-2layer"
figdir = home * "/figures"
datadir = home * "/data"

# for plotting
latex_fonts = Theme(major_label_font="CMU Serif", major_label_font_size=16pt,
                    minor_label_font="CMU Serif", minor_label_font_size=14pt,
                    key_title_font="CMU Serif", key_title_font_size=12pt,
                    key_label_font="CMU Serif", key_label_font_size=10pt);
gadfly_colors = Scale.color_discrete_hue().f(10)
my_color_subset = [1 6 4 5 2 7];
my_colors = gadfly_colors[my_color_subset];
my_colormap = p -> get(ColorSchemes.linear_blue_5_95_c73_n256,p)
scatter_theme = Theme(discrete_highlight_color = c -> nothing,alphas = [.5],point_size = 1pt,line_width = 0pt);``

# TYPES

mutable struct AdvectionDiffusion
    K1::Function
    dK1dx::Function
    U1::Function
    reflecting::Bool
end

# mutable struct PathFunctions
#     δ::Function
#     g::Function
# end

mutable struct PathFunctions2
    δ::Function
    g::Function
    Ky::Function
end

# mutable struct ItoProcess
#     a::Function
#     b::Function
#     K1::Function
#     dK1dx::Function
#     U1::Function
#     reflecting::Bool
#     δ::Function
#     g::Function
# end

mutable struct ItoProcess2
    a::Function
    b::Function
    K1::Function
    dK1dx::Function
    U1::Function
    reflecting::Bool
    δ::Function
    g::Function
    Ky::Function
end

# struct LarvaOutput
#     sim::DataFrame
#     stats::DataFrame
#     time::Float64
# end

mutable struct Larva2Output
    sim::DataFrame
    stats::DataFrame
    F::ItoProcess2
    time::Float64
end


# struct LarvaeOutput
#     sims::DataFrame
#     stats::DataFrame
#     stats_settled::DataFrame
#     meta::DataFrame
#     time::Float64
# end

struct Larvae2Output
    sims::DataFrame
    stats::DataFrame
    stats_settled::DataFrame
    meta::DataFrame
    F::ItoProcess2
    time::Float64
end

mutable struct VSwim
    vm::Function
    stochastic::Bool
    args::Vector # besides p
end


# struct LarvaSweepOutput
#     sims::DataFrame
#     scores::DataFrame
#     scores_settled::DataFrame
#     meta::DataFrame
#     flow_dict::OrderedDict
#     vm_dict::OrderedDict
#     pars::OrderedDict
#     parameter_combos::DataFrame
#     runTime::Float64
# end

# SUPPORTING FUNCTIONS

function apply_to_grid(f::Function,x1::Array{Float64,1},x2::Array{Float64,1}; names = (:x1,:x2,:f),as_dataframe = true)
    N1 = length(x1)
    N2 = length(x2)
    F = zeros(N2,N1)
    for n1 = 1:N1
        for n2 = 1:N2
            F[n2,n1] = f(x1[n1],x2[n2])
        end
    end
    if as_dataframe
        name1 = names[1]
        name2 = names[2]
        name3 = names[3]
        F = DataFrame(
            name1 => repeat(x1,inner = N2),
            name2 => repeat(x2,outer = N1),
            name3 => reshape(F,N1*N2)
        )
    end
    return(F)
end

function create_contour_layers(df::DataFrame,z::Symbol,c::Float64;color = "white",line_width = 1mm)
    # for now this is specific to this application :( no choice of x,y
    # Avec = unique(df.A); LA = length(Avec)
    # Bvec = unique(df.B); LB = length(Bvec)
    Avec = unique(df.A); LA = length(Avec)
    Bvec = unique(df.B); LB = length(Bvec)
    Z = zeros(LA,LB)
    for i in 1:LA
        for j in 1:LB
            ℓ = LB*(i-1)+j
            Z[i,j] = df[ℓ,z]
        end
    end
    CC = lines(contour(Avec,Bvec,Z,c))
    coords = [coordinates(C) for C in CC]
    layers = [layer(x = xy[1],y = xy[2],Geom.line(preserve_order = true),Theme(default_color = color,line_width = line_width)) for xy in coords]
    return(layers)
end

function initialize_df(Avec,Bvec,cols)
    df = DataFrame(:A => repeat(Avec,inner = NB),:B => repeat(Bvec,outer = NA));
    for col in cols
        df[!,col] .= 0.
    end
    return(df)
end



# MODEL FUNCTIONS
function EulerMaruyama(p,t,x,z,F,db)
    return(x + F.a(p,t,x,z)*p[:dt] + F.b(p,t,x,z)*db)
end

function RungeKutta(p,t,x,z,F,db)
    dt_sq = sqrt(p[:dt])
    χ = x + F.a(p,t,x,z)*p[:dt] + F.b(p,t,x,z)*dt_sq
    return(x + F.a(p,t,x,z)*p[:dt] + F.b(p,t,x,z)*db + .5(F.b(p,t,χ,z) - F.b(p,t,x,z))*(db^2 - p[:dt])/dt_sq)
end

include("parameters2.jl")
include("fx2.jl")
include("Larva3.jl")
include("Larvae3.jl")
