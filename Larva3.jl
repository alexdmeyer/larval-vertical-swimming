
function Larva(
    p::Dict,
    F::AdvectionDiffusion,
    pathfx::PathFunctions2,
    V::Array{Int64,1};
    sim::DataFrame = DataFrame(
        :t => p[:tvec],
        :X => p[:X0]*ones(p[:NT]),
        :Z => V,
        :E_plus => zeros(p[:NT]),
        :E_minus => zeros(p[:NT]),
        :D => zeros(p[:NT])
        ),
    no_settling::Bool = false,
    dB::Array{Float64,1} = randn(p[:NT]-1)*sqrt(p[:dt]),
    solver::Function = EulerMaruyama
    )


    t0 = time()

    # convert F to an Ito, with all the fixins
    F = AD_to_Ito(F,pathfx)

    dt = p[:dt]
    dt_sq = sqrt(p[:dt])

    flag = 0 # becomes 1 if we settle
    n_final = p[:NT] # index of last time

    for n = 1:(p[:NT]-1)
        x = sim.X[n]
        z = sim.Z[n]
        Δz = sim.Z[n+1]-z
        t = sim.t[n]

        xx = solver(p,t,x,z,F,dB[n])
        if xx >= 0.
            sim.X[n+1] = xx
        else
            xx = ifelse(F.reflecting,-xx,0.)
            sim.X[n+1] = xx
        end

        sim.E_plus[n+1] = sim.E_plus[n] + F.g(p,t,x,z)*dt
        sim.E_minus[n+1] = sim.E_minus[n] + p[:m]*dt + p[:C][Δz]
        sim.D[n+1] = sim.D[n] + F.δ(p,t,x,z)*dt

        # did we settle?
        if t+dt >= p[:T_PC] && xx <= p[:xh]
            flag = 1
            if !no_settling
                n_final = n+1
                break
            end
        end
    end

    stats = DataFrame(
        :settled => flag,
        :risk => sim.D[n_final],
        :p_surv => exp(-sim.D[n_final]),
        :E_gain => sim.E_plus[n_final],
        :E_loss => sim.E_minus[n_final],
        :XT => sim.X[n_final],
        :T => p[:tvec][n_final],
        :alongshore_movement => sqrt(2*sum([pathfx.Ky(p,sim.X[i],sim.Z[i]) for i in 1:n_final])*p[:dt]),
        :k => n_final
    )

    return(Larva2Output(sim,stats,F,time() - t0))
end


#
# function inner_sim(
#     p::Dict,
#     F::ItoProcess2,
#     V::Array{Int64,1},
#     solver::Function,
#     dB::Array{Float64,1};
#     no_settling::Bool = false,
#     sim::DataFrame = DataFrame(:t => p[:tvec],:X => p[:X0]*ones(p[:NT]),:Z => V,:E_plus => zeros(p[:NT]),:E_minus => zeros(p[:NT]),:D => zeros(p[:NT]))
#     )
#
#     dt = p[:dt]
#     dt_sq = sqrt(p[:dt])
#
#     flag = 0 # becomes 1 if we settle
#     n_final = p[:NT] # index of last time
#
#     for n = 1:(p[:NT]-1)
#         x = sim.X[n]
#         z = sim.Z[n]
#         Δz = sim.Z[n+1]-z
#         t = sim.t[n]
#
#         xx = solver(p,t,x,z,F,dB[n])
#         # χ = x + F.a(p,t,x,z)*dt + F.b(p,t,x,z)*dt_sq
#         # xx = x + F.a(p,t,x,z)*dt + F.b(p,t,x,z)*dB[n] + .5(F.b(p,t,χ,z) - F.b(p,t,x,z))*(dB[n]^2 - dt)/dt_sq
#         if xx >= 0.
#             sim.X[n+1] = xx
#         else
#             xx = ifelse(F.reflecting,-xx,0.)
#             sim.X[n+1] = xx
#         end
#
#         sim.E_plus[n+1] = sim.E_plus[n] + F.g(p,t,x,z)*dt
#         sim.E_minus[n+1] = sim.E_minus[n] + p[:m]*dt + p[:C][Δz]
#         sim.D[n+1] = sim.D[n] + F.δ(p,t,x,z)*dt
#
#         # did we settle?
#         if t+dt >= p[:T_PC] && xx <= p[:xh]
#             flag = 1
#             if !no_settling
#                 n_final = n+1
#                 break
#             end
#         end
#     end
#
#     return(sim,flag,n_final)
# end
#
#
