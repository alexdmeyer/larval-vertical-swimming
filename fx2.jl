# TYPES

function AD_to_Ito(F::AdvectionDiffusion,pathfx::PathFunctions2)

    function a(p,t,x,z)
        return(F.U1(p,t,x,z) + F.dK1dx(p,t,x,z))
    end

    function b(p,t,x,z)
        return(sqrt(2*F.K1(p,t,x,z)))
    end

    return(ItoProcess2(a,b,F.K1,F.dK1dx,F.U1,F.reflecting,pathfx.δ,pathfx.g,pathfx.Ky))
end

# the flow scheme

function K(p::Dict,_,_,z) # p,t,x,z
    return(ifelse(z == 1,p[:Kx1],p[:Kx0]))
    # return(p[:k_inf] * y[1] / (p[:x_half] + y[1]))
end;

function ∂x_K(p::Dict,_,_,_)
    return(0.)
    # return(p[:k_inf] * p[:x_half] / (p[:x_half] + y[1])^2.)
end;

function U(p::Dict,_,_,z)
    return(ifelse(z == 0.,p[:U_minus],p[:U_plus]))
end

function K_cbl(p::Dict,_,x,_)
    # return(p[:K0])
    return(p[:K0] * x / (p[:x_half] + x))
end;

function ∂x_K_cbl(p::Dict,_,x,_)
    # return(0.)
    return(p[:K0] * p[:x_half] / (p[:x_half] + x)^2.)
end;

function U_cbl(p::Dict,_,x,z)
    return(ifelse(z == 0,p[:U_minus],p[:U_plus]) * x / (p[:x_half] + x))
end

UW = AdvectionDiffusion(K,∂x_K,U,true)
UW_cbl = AdvectionDiffusion(K_cbl,∂x_K_cbl,U_cbl,false)

function vm_generic2(p,A,B)
    V = zeros(p[:NT])
    V[.!(A/2 .<= (p[:tvec] .% 24.) .< (24 - A/2)) .& (p[:tvec] .< B*p[:T_LD])] .= 1
    return(Int.(V))
end

function vm_generic2(V,p,A,B)
    V[.!(A/2 .<= (p[:tvec] .% 24.) .< (24 - A/2)) .& (p[:tvec] .< B*p[:T_LD])] .= 1
    return(Int.(V))
end

function vm_neutral2(p,λ1,λ2) # \lambda = mean length of stay in lower layer
    V = zeros(p[:NT])
    ℓ = 0
    t = 0.
    i = 1
    while t < p[:T_LD]
        i += 1
        tnext = t + ifelse(iseven(i),λ1,λ2)*randexp()
        V[t .<= p[:tvec] .< tnext] .= ℓ
        ℓ = 1 - ℓ
        t = tnext
    end
    return(Int.(V))
end

function vm_neutral2(V,p,λ1,λ2) # \lambda = mean length of stay in lower layer
    ℓ = 0
    t = 0.
    i = 1
    while t < p[:T_LD]
        i += 1
        tnext = t + ifelse(iseven(i),λ1,λ2)*randexp()
        V[t .<= p[:tvec] .< tnext] .= ℓ
        ℓ = 1 - ℓ
        t = tnext
    end
    return(Int.(V))
end

# swims = OrderedDict(
#     :Neutral => VSwim(vm_neutral2,true,[18.,6.]),
#     :Binary => VSwim(vm_generic2,false,[24.,.2]),
#     :DVM => VSwim(vm_generic2,false,[6.,1.]),
#     :Hybrid => VSwim(vm_generic2,false,[6,.5])
# );

swims = OrderedDict(
    :PASSIVE => VSwim(vm_neutral2,true,[14.,1.]),
    :BINARY => VSwim(vm_generic2,false,[24.,.2]),
    :DVM1 => VSwim(vm_generic2,false,[6.,1.]),
    :DVM2 => VSwim(vm_generic2,false,[9.,1.]),
    :HYBRID => VSwim(vm_generic2,false,[6,.5])
);

# # PATH FUNCTIONS

pathfx2 = PathFunctions2(
    (p,t,x,z) -> ifelse(z == 0 || !(p[:h_light] <= (t % 24.) < p[:h_dark]),p[:δ0],p[:δ1]),
    (p,t,x,z) -> ifelse(z == 0,p[:g0],p[:g1]),
    (p,x,z) -> ifelse(z == 0,p[:Ky0],p[:Ky1]) * x / (p[:ky] + x)
)

pathfx2h = PathFunctions2(
    (p,t,x,z) -> ifelse(x <= p[:xδ],p[:δ1],p[:δ0]),
    (p,t,x,z) -> ifelse((p[:xg2] > x > p[:xg1]) && (z == 1.),p[:g1],p[:g0]),
    (p,x,z) -> ifelse(z == 0,p[:Ky0],p[:Ky1]) * x / (p[:ky] + x)
)
