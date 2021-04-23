function Larvae(
    p::Dict,
    F::AdvectionDiffusion,
    pathfx::PathFunctions2,
    VS::VSwim,
    NSIM::Int64;
    lite::Bool = true,
    no_settling::Bool = false,
    solver::Function = EulerMaruyama
    )

    t0 = time()

    # draw random numbers
    DB = randn(p[:NT]-1,NSIM)*sqrt(p[:dt])

    # the first sim to initialize DF
    V = VS.vm(p,VS.args...)
    L = Larva(p,F,pathfx,V,no_settling = no_settling,solver = solver)

    stats = DataFrame(eltypes(L.stats),names(L.stats),NSIM)
    stats[1,:] = L.stats

    if lite
        sims = DataFrame()
    else
        sims = DataFrame(eltypes(L.sim),names(L.sim),p[:NT]*NSIM)
        sims[1:p[:NT],:] = L.sim
    end

    Threads.@threads for n = 2:NSIM
        if VS.stochastic
            V = VS.vm(V,p,VS.args...)
            # L.sim.Z = V
        end
        L = Larva(p,F,pathfx,V,dB = DB[:,n],no_settling = no_settling,solver = solver) # ,sim = L.sim)
        stats[n,:] = L.stats
        if !lite
            sims[(n-1)*p[:NT] .+ (1:p[:NT]),:] = L.sim
        end
    end

    stats[!,:n] = 1:NSIM
    stats_settled = filter(r -> r.settled == 1,stats)
    sims[!,:n] = repeat(1:NSIM,inner = p[:NT])


    meta = DataFrame(
        :frac_settled => sum(stats.settled)/NSIM,
        :mean_risk => mean(stats_settled.risk),
        :mean_p_surv => mean(stats_settled.p_surv), # meaningless for non-settled larvae
        :mean_E_gain => mean(stats_settled.E_gain),
        :mean_E_loss => mean(stats_settled.E_loss),
        :mean_alongshore_movement => mean(stats_settled.alongshore_movement),
        :mean_settle_time => mean(stats_settled.T)
    )
    return(Larvae2Output(sims,stats,stats_settled,meta,L.F,time() - t0))
end


function Rescore(
    p::Dict,
    L::Larvae2Output,
    pathfx::PathFunctions2,
    )

    if nrow(L.sims) == 0
        error("I need some simulations!")
    end
    if (p[:dt] != L.sims.t[2] - L.sims.t[1]) || (p[:T_LD] != maximum(L.sims.t))
        error("You can't rescore if you've changed your time vector!")
    end

    NSIM = L.sims.n[end]
    sims = copy(L.sims)
    stats = copy(L.stats)
    meta = copy(L.meta)

    # Threads.@threads for n = 1:NSIM
    for n = 1:NSIM
        X = L.sims[(n-1)*p[:NT] .+ (1:p[:NT]),:X]
        Z = L.sims[(n-1)*p[:NT] .+ (1:p[:NT]),:Z]
        k = L.stats.k[n]
        K = k + (n-1)*p[:NT]
        for M = 1:(k-1)
            m = M + (n-1)*p[:NT]
            x = X[M]
            z = Z[M]
            t = p[:tvec][M]
            sims[m+1,:E_plus] = sims[m,:E_plus] + pathfx.g(p,t,x,z)*p[:dt]
            sims[m+1,:E_minus] = sims[m,:E_minus] + p[:m]*p[:dt] + p[:C][Z[M+1]-z]
            sims[m+1,:E_net] = sims[m+1,:E_plus] - sims[m+1,:E_minus]
            sims[m+1,:D] = sims[m,:D] + pathfx.δ(p,t,x,z) * p[:dt]
        end

        stats.risk[n] = sims.D[K]
        stats.p_surv[n] = exp(-sims.D[K])
        stats.E_gain[n] = sims.E_plus[K]
        stats.E_loss[n] = sims.E_minus[K]
        stats.E_ratio[n] = sims.E_plus[K]/sims.E_minus[K]
        stats.ΔE[n] = sims.E_net[K]
        stats.alongshore_movement[n] = sum([pathfx.Ky(p,X[i],Z[i]) for i in 1:k])*p[:dt]

    end

    stats_settled = filter(r -> r.settled == 1,stats)
    F = ItoProcess2(L.F.a,L.F.b,L.F.K1,L.F.dK1dx,L.F.U1,L.F.reflecting,pathfx.δ,pathfx.g,pathfx.Ky)
    meta[1,:mean_p_surv] = mean(stats_settled.p_surv)
    meta[1,:mean_E_gain] = mean(stats_settled.E_gain)
    meta[1,:mean_E_loss] = mean(stats_settled.E_loss)
    meta[1,:mean_alongshore_movement] = mean(stats_settled.alongshore_movement)
    meta[1,:mean_risk] = mean(stats_settled.risk)


    return(Larvae2Output(sims,stats,stats_settled,meta,F,L.time))
end

function Rescore_energy(p::Dict,df::DataFrame)
    X = copy(df)
    N = nrow(X)
    for n = 1:N
        A = X.A[n]
        B = X.B[n]
        T = X.mean_settle_time[n]
        if A < 24.
            N_down = floor(B*T/24) + ifelse((B*T % 24) > .5*A,1,0)
            N_up = floor(B*T/24) + ifelse((B*T % 24) > 24 - .5*A,1,0)
        else
            N_down = ifelse(B < 1,1,0)
            N_up = 0
        end
        X.mean_E_loss[n] = p[:m]*T + p[:C][-1]*N_down + p[:C][1]*N_up
    end
    return(X)
end
