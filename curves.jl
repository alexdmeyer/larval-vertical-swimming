include("requirements2.jl")
NSIM2 = 2000
p = parameters2()
swims = OrderedDict(
    :Neutral => VSwim(vm_neutral2,true,[14.,1.]),
    :Binary => VSwim(vm_generic2,false,[24.,.25]),
    :DVM => VSwim(vm_generic2,false,[10.,1.]),
    :Hybrid => VSwim(vm_generic2,false,[10,.25])
);
Nswim = length(swims)

Tvec = 24. * [5,25,50,100]; NLD = length(Tvec);
Kvec = 10. .^ [-2,-1,0,1,2]; NK = length(Kvec);
Uplusvec = collect(0:.2:1); NU = length(Uplusvec);

ex_T = DataFrame(
    :T_LD => repeat(Tvec,outer = Nswim),
    :swim => repeat(collect(keys(swims)),inner = NLD),
    :frac_settled => zeros(Nswim*NLD),
    :mean_alongshore_movement => zeros(Nswim*NLD),
    :mean_settle_time => zeros(Nswim*NLD),
);

t0 = time();
for j = 1:(NLD*Nswim)
    this_p = parameters2(T_LD = ex_T.T_LD[j])
    L = Larvae(this_p,UW,pathfx2h,swims[ex_T.swim[j]],NSIM2)
    ex_T.frac_settled[j] = L.meta.frac_settled[1]
    ex_T.mean_alongshore_movement[j] = L.meta.mean_alongshore_movement[1]
    ex_T.mean_settle_time[j] = L.meta.mean_settle_time[1]
end
t1 = time() - t0

ex_K = DataFrame(
    :K => repeat(Kvec,outer = Nswim),
    :swim => repeat(collect(keys(swims)),inner = NK),
    :frac_settled => zeros(Nswim*NK),
    :mean_alongshore_movement => zeros(Nswim*NK)
);
t0 = time();
for j = 1:(NK*Nswim)
    this_p = parameters2(T_LD = 24*25.);
    this_p[:Kx1] = ex_K.K[j]
    this_p[:Kx0] = .25*ex_K.K[j]
    L = Larvae(this_p,UW,pathfx2,swims[ex_K.swim[j]],NSIM2)
    ex_K.frac_settled[j] = L.meta.frac_settled[1]
    ex_K.mean_alongshore_movement[j] = L.meta.mean_alongshore_movement[1]
end
t1 = time() - t0

ex_U = DataFrame(
    :Uplus => repeat(Uplusvec,outer = Nswim),
    :swim => repeat(collect(keys(swims)),inner = NU),
    :frac_settled => zeros(Nswim*NU),
);

t0 = time();
for j = 1:(NU*Nswim)
    this_p = parameters2(T_LD = 24*25.);
    this_p[:U_plus] = ex_U.Uplus[j]
    this_p[:U_minus] = -.25*ex_U.Uplus[j]
    L = Larvae(this_p,UW,pathfx2,swims[ex_U.swim[j]],NSIM2)
    ex_U.frac_settled[j] = L.meta.frac_settled[1]
end
t1 = time() - t0

d1 = plot(ex_T,x = :T_LD,y = :frac_settled,color = :swim,Geom.point,Geom.line,
    Guide.xlabel("Larval duration, T (h)"),
    Guide.ylabel("Frac settled"),
    Guide.colorkey("Behavior"),
    latex_fonts
)

e1 = plot(ex_U,x = :Uplus,y = :frac_settled,color = :swim,Geom.point,Geom.line,
    Guide.xlabel("Velocity of offshore layer (km/h)"),
    Guide.ylabel("Frac settled"),
    Guide.colorkey("Behavior"),
    latex_fonts
)

f1 = plot(ex_K,x = :K,y = :frac_settled,color = :swim,Geom.point,Geom.line,
    Scale.x_log10,
    Guide.xlabel("Diffusivity (km^2/h)"),
    Guide.ylabel("Frac settled"),
    Guide.colorkey("Behavior"),
    latex_fonts,
)


# draw(PDF(figdir * "/fig1/frac_settled_short.pdf",13cm,10cm),a1)
# draw(PDF(figdir * "/fig1/frac_settled_long.pdf",13cm,10cm),b1)
# draw(PDF(figdir * "/fig1/mean_settle_time.pdf",13cm,10cm),c1)
draw(PDF(figdir * "/sensitivity/frac_settled_v_T.pdf",13cm,10cm),d1)
draw(PDF(figdir * "/sensitivity/frac_settled_v_U.pdf",13cm,10cm),e1)
draw(PDF(figdir * "/sensitivity/frac_settled_v_K.pdf",13cm,10cm),f1)
