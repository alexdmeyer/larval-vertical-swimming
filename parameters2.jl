function parameters2(;T_LD::Float64 = 24. * 25.,pc::Float64 = 0.6,K0::Float64 = 50 * 3.6e-3, layer_ratio::Float64 = 4.,U_minus::Float64 = -.02*3.6)
    p = Dict()

    dt = 1//4
    p[:dt] = dt      # time step

    # larval duration
    p[:T_LD] = T_LD
    p[:T_PC] = pc * T_LD # pre-competence duration
    p[:tvec] = collect(0:dt:T_LD)
    p[:NT] = length(0:dt:T_LD)

    # physical structure
    p[:X0] = 0.5  # initial X
    p[:xh] = 5. # habitat width
    p[:layer_ratio] = layer_ratio # how much thicker is top layer than bottom layer?

    # light/dark cycle
    p[:h_dark] = 18.    # dark at 6pm
    p[:h_light] = 6.    # light at 6am

    # horizontal flow
    p[:x_half] = 10.   # offshore dist at which half max speed attained
    p[:Kx0] = K0 # 200*3.6e-3
    p[:Kx1] = layer_ratio*K0 # 200*3.6e-3
    p[:U_minus] = U_minus
    p[:U_plus] = -layer_ratio*U_minus
    p[:ky] = 10.
    p[:Ky0] = K0
    p[:Ky1] = layer_ratio*K0


    # energy
    p[:m] = 15 # maintenance energy rate
    p[:C] = Dict(-1 => 0.,0 => 0.,1 => 10.) # cost of switch: 0 = low -> high, 1 = high -> low
    p[:g0] = 1e4          # phytoplankton concentration offshore + at depth
    p[:g1] = 1e5 # conx in upper layer
    # p[:gc] = 1e3         # phytoplankton conx at coast + at depth
    # p[:gs] = 1e5          # phytoplankton conx at surface + offshore
    # p[:f] = 2.5e-6          # rate of energy uptake per unit phytoplankton conx
    p[:xg1] = 25.
    p[:xg2] = 75.

    # mortality risk
    p[:δ0] = .005          # mortality rate during darkness
    p[:δ1] = .025           # mortality rate at surface + during light
    p[:xδ] = 5.          # dist at which coast effect is halved
    # p[:δs] = .025          # mortality rate at surface during daylight + offshore

    return(p)
end
