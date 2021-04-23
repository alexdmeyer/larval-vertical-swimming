include("requirements2.jl")
NSIM = 5000
p = parameters2()
cols = [:frac_settled,:mean_alongshore_movement,:mean_p_surv,:mean_risk,:mean_E_gain,:mean_E_loss,:mean_settle_time]

p0 = parameters2(T_LD = 120.)
pT = parameters2(T_LD = 360.) # slightly smaller
pU = parameters2(U_minus = -.03*3.6) # slightly bigger
pK = parameters2(K0 = 75*3.6e-3) # slightly bigger

LA = 49
LB = 21
Avec = range(0.,24.,length = LA)
Bvec = range(0.,1.,length = LB)
frontier = DataFrame(:A => Avec[2:end],:B0 => .5*(p[:U_plus]/p[:U_minus] .+ sqrt.((p[:U_plus]/p[:U_minus])^2 .+ 96 ./ Avec[2:end])))
examples = DataFrame(:behavior => [:OVM,:DVM,:HYBRID],:A => [24,10,10],:B => [.25,1,.25])

v = CSV.read(datadir*"/v-600.csv")
h = CSV.read(datadir*"/h-600.csv")
v0 = CSV.read(datadir*"/v-120.csv")
vT = CSV.read(datadir*"/v-360.csv")
vU = CSV.read(datadir*"/v-600U.csv")
vK = CSV.read(datadir*"/v-600K.csv")

# v = DataFrame(:A => repeat(Avec,inner = LB),:B => repeat(Bvec,outer = LA));
# for col in cols
#     v[!,col] .= 0.
# end
# h = copy(v);
# v0 = copy(v);
# vT = copy(v);
# vU = copy(v);
# vK = copy(v);
#
# k = 0
# Threads.@threads for n = 1:LA*LB
#
#     VS = VSwim(vm_generic2,false,[v.A[n],v.B[n]])
#
#     # the standard comp, vertical and horizontal
#     L = Larvae(p,UW,pathfx2,VS,NSIM,lite = true)
#     v[n,cols] = L.meta[1,cols]
#     L = Larvae(p,UW,pathfx2h,VS,NSIM,lite = true)
#     h[n,cols] = L.meta[1,cols]
#
#     # short T
#     L = Larvae(p0,UW,pathfx2,VS,NSIM,lite = true)
#     v0[n,cols] = L.meta[1,cols]
#
#     # perturbed T
#     L = Larvae(pT,UW,pathfx2,VS,NSIM,lite = true)
#     vT[n,cols] = L.meta[1,cols]
#
#     # perturbed U
#     L = Larvae(pU,UW,pathfx2,VS,NSIM,lite = true)
#     vU[n,cols] = L.meta[1,cols]
#
#     # perturbed K
#     L = Larvae(pK,UW,pathfx2,VS,NSIM,lite = true)
#     vK[n,cols] = L.meta[1,cols]
#
#     global k += 1
#     println(string(100k/(LA*LB)) * " percent complete")
#
# end
# #
#
v[!,:EM] = p[:m]*v.mean_settle_time
v[!,:EL] = v.mean_E_loss - v.EM

v[!,:sensitivity_T] = ifelse.(v.frac_settled .> 0,sign.(v.frac_settled - vT.frac_settled),NaN)
v[!,:sensitivity_U] = ifelse.(v.frac_settled .> 0,sign.(-v.frac_settled + vU.frac_settled),NaN)
v[!,:sensitivity_K] = ifelse.(v.frac_settled .> 0,sign.(-v.frac_settled + vK.frac_settled),NaN)

V = filter(r -> r.frac_settled > 0,v)
H = filter(r -> r.frac_settled > 0,h)
#
# # compute contour values for neutral behaviors
# # add contours for Neutral behavior....
# PassiveFloat = VSwim(vm_neutral2,true,[14.,1.])
# c0 = Larvae(p0,UW,pathfx2,PassiveFloat,NSIM,lite = false);
# c0.meta[1,:mean_E_loss] = p0[:m]*c0.meta.mean_settle_time[1];
# c = Larvae(p,UW,pathfx2,PassiveFloat,NSIM,lite = false);
# c.meta[1,:mean_E_loss] = p[:m]*c.meta.mean_settle_time[1];
# ch = Larvae(p,UW,pathfx2h,PassiveFloat,NSIM,lite = false);
# ch.meta[1,:mean_E_loss] = p[:m]*ch.meta.mean_settle_time[1];
#
#
#
# CSV.write(datadir * "/v-600.csv",v)
# CSV.write(datadir * "/h-600.csv",h)
# CSV.write(datadir * "/v-120.csv",v0)
# CSV.write(datadir * "/v-360.csv",vT)
# CSV.write(datadir * "/v-600U.csv",vU)
# CSV.write(datadir * "/v-600K.csv",vK)

VV = copy(V);
VV[!,:A2] = VV.A .+ (Avec[2] - Avec[1])
VV[!,:B2] = VV.B .+ (Bvec[2] - Bvec[1])

sensU = plot(
    layer(examples,x = :A,y = :B,shape = :behavior,Geom.point,Theme(default_color = colorant"black",point_size = 1.5mm)),
    layer(frontier,x = :A,y = :B0,Geom.line,Theme(line_style = [:dash],default_color = colorant"gray",line_width=1mm)),
    layer(filter(r -> r.sensitivity_U == -1,VV),xmin = :A,ymin = :B,xmax = :A2,ymax = :B2,color = [my_colors[1]],Geom.rect),
    layer(filter(r -> r.sensitivity_U == 1,VV),xmin = :A,ymin = :B,xmax = :A2,ymax = :B2,color = [my_colors[2]],Geom.rect),
    layer(filter(r -> r.sensitivity_U == 0,VV),xmin = :A,ymin = :B,xmax = :A2,ymax = :B2,color = [my_colors[6]],Geom.rect),
    # layer(V,x = :A,y = :B,color = :sensitivity_T,Geom.rectbin),
    Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Guide.xlabel("A (hours at surface)"),
    Guide.ylabel("B (timing of switch)"),
    # Guide.title(""),
    Guide.colorkey(""),
    latex_fonts
)

sensT = plot(
    layer(examples,x = :A,y = :B,shape = :behavior,Geom.point,Theme(default_color = colorant"black",point_size = 1.5mm)),
    layer(frontier,x = :A,y = :B0,Geom.line,Theme(line_style = [:dash],default_color = colorant"gray",line_width=1mm)),
    layer(filter(r -> r.sensitivity_T == -1,VV),xmin = :A,ymin = :B,xmax = :A2,ymax = :B2,color = [my_colors[1]],Geom.rect),
    layer(filter(r -> r.sensitivity_T == 1,VV),xmin = :A,ymin = :B,xmax = :A2,ymax = :B2,color = [my_colors[2]],Geom.rect),
    layer(filter(r -> r.sensitivity_T == 0,VV),xmin = :A,ymin = :B,xmax = :A2,ymax = :B2,color = [my_colors[6]],Geom.rect),
    # layer(V,x = :A,y = :B,color = :sensitivity_T,Geom.rectbin),
    Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Guide.xlabel("A (hours at surface)"),
    Guide.ylabel("B (timing of switch)"),
    # Guide.title(""),
    Guide.colorkey(""),
    latex_fonts
)

sensK = plot(
    layer(examples,x = :A,y = :B,shape = :behavior,Geom.point,Theme(default_color = colorant"black",point_size = 1.5mm)),
    layer(frontier,x = :A,y = :B0,Geom.line,Theme(line_style = [:dash],default_color = colorant"gray",line_width=1mm)),
    layer(filter(r -> r.sensitivity_K == -1,VV),xmin = :A,ymin = :B,xmax = :A2,ymax = :B2,color = [my_colors[1]],Geom.rect),
    layer(filter(r -> r.sensitivity_K == 1,VV),xmin = :A,ymin = :B,xmax = :A2,ymax = :B2,color = [my_colors[2]],Geom.rect),
    layer(filter(r -> r.sensitivity_K == 0,VV),xmin = :A,ymin = :B,xmax = :A2,ymax = :B2,color = [my_colors[6]],Geom.rect),
    # layer(V,x = :A,y = :B,color = :sensitivity_T,Geom.rectbin),
    Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Guide.xlabel("A (hours at surface)"),
    Guide.ylabel("B (timing of switch)"),
    # Guide.title(""),
    Guide.colorkey(""),
    latex_fonts
)

draw(PDF(figdir*"/sensitivity/heatmap-U.pdf",13cm,10cm),sensU)
draw(PDF(figdir*"/sensitivity/heatmap-K.pdf",13cm,10cm),sensK)
draw(PDF(figdir*"/sensitivity/heatmap-T.pdf",13cm,10cm),sensT)


DE0 = plot(
    create_contour_layers(v0,:frac_settled,c0.meta[1,:frac_settled])...,
    layer(examples,x = :A,y = :B,Geom.point,shape = :behavior,Theme(default_color = colorant"black",point_size = 1.5mm)),
    layer(frontier,x = :A,y = :B0,Geom.line,Theme(line_style = [:dash],default_color = colorant"gray",line_width=1mm)),
    layer(v0,x = :A,y = :B,color = :frac_settled,Geom.rectbin),
    Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = p -> get(ColorSchemes.linear_blue_5_95_c73_n256,p)),
    Guide.xlabel("A (hours at surface)"),
    Guide.ylabel("B (timing of switch)"),
    # Guide.title("Dispersal Efficiency, DE (T = 5d)"),
    Guide.colorkey(""),
    latex_fonts
)

DE = plot(
    create_contour_layers(v,:frac_settled,c.meta[1,:frac_settled])...,
    layer(examples,x = :A,y = :B,Geom.point,shape = :behavior,Theme(default_color = colorant"black",point_size = 1.5mm)),
    layer(frontier,x = :A,y = :B0,Geom.line,Theme(line_style = [:dash],default_color = colorant"gray",line_width=1mm)),
    layer(v,x = :A,y = :B,color = :frac_settled,Geom.rectbin),
    Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = p -> get(ColorSchemes.linear_blue_5_95_c73_n256,p)),
    Guide.xlabel("A (hours at surface)"),
    Guide.ylabel("B (timing of switch)"),
    # Guide.title("Dispersal Efficiency, DE (T = 25d)"),
    Guide.colorkey(""),
    latex_fonts
)


DT = plot(
    create_contour_layers(v,:mean_settle_time,c.meta[1,:mean_settle_time])...,
    layer(examples,x = :A,y = :B,Geom.point,shape = :behavior,Theme(default_color = colorant"black",point_size = 1.5mm)),
    layer(frontier,x = :A,y = :B0,Geom.line,Theme(line_style = [:dash],default_color = colorant"gray",line_width=1mm)),
    layer(V,x = :A,y = :B,color = :mean_settle_time,Geom.rectbin),
    Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = p -> get(ColorSchemes.linear_blue_5_95_c73_n256,p)),
    Guide.xlabel("A (hours at surface)"),
    Guide.ylabel("B (timing of switch)"),
    # Guide.title("Mean Dispersal Time, ET*"),
    Guide.colorkey(""),
    latex_fonts
)

DP = plot(
    create_contour_layers(v,:mean_alongshore_movement,c.meta[1,:mean_alongshore_movement])...,
    layer(examples,x = :A,y = :B,Geom.point,shape = :behavior,Theme(default_color = colorant"black",point_size = 1.5mm)),
    layer(frontier,x = :A,y = :B0,Geom.line,Theme(line_style = [:dash],default_color = colorant"gray",line_width=1mm)),
    layer(V,x = :A,y = :B,color = :mean_alongshore_movement,Geom.rectbin),
    Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = p -> get(ColorSchemes.linear_blue_5_95_c73_n256,p)),
    Guide.xlabel("A (hours at surface)"),
    Guide.ylabel("B (timing of switch)"),
    # Guide.title("Dispersal Potential, DP"),
    Guide.colorkey(""),
    latex_fonts
)

S_diurnal = plot(
    create_contour_layers(v,:mean_p_surv,c.meta[1,:mean_p_surv])...,
    layer(examples,x = :A,y = :B,Geom.point,shape = :behavior,Theme(default_color = colorant"black",point_size = 1.5mm)),
    layer(frontier,x = :A,y = :B0,Geom.line,Theme(line_style = [:dash],default_color = colorant"gray",line_width=1mm)),
    layer(V,x = :A,y = :B,color = :mean_p_surv,Geom.rectbin),
    Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = p -> get(ColorSchemes.linear_blue_5_95_c73_n256,p)),
    Guide.xlabel("A (hours at surface)"),
    Guide.ylabel("B (timing of switch)"),
    # Guide.title("Prob of survival, S (diurnal)"),
    Guide.colorkey(""),
    latex_fonts
)

S_nearshore = plot(
    create_contour_layers(h,:mean_p_surv,ch.meta[1,:mean_p_surv])...,
    layer(examples,x = :A,y = :B,Geom.point,shape = :behavior,Theme(default_color = colorant"black",point_size = 1.5mm)),
    layer(frontier,x = :A,y = :B0,Geom.line,Theme(line_style = [:dash],default_color = colorant"gray",line_width=1mm)),
    layer(H,x = :A,y = :B,color = :mean_p_surv,Geom.rectbin),
    Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = p -> get(ColorSchemes.linear_blue_5_95_c73_n256,p),maxvalue = .1),
    Guide.xlabel("A (hours at surface)"),
    Guide.ylabel("B (timing of switch)"),
    # Guide.title("Prob of survival, S (nearshore)"),
    Guide.colorkey(""),
    latex_fonts
)

EL = plot(
    layer(examples,x = :A,y = :B,Geom.point,shape = :behavior,Theme(default_color = colorant"black",point_size = 1.5mm)),
    layer(frontier,x = :A,y = :B0,Geom.line,Theme(line_style = [:dash],default_color = colorant"gray",line_width=1mm)),
    layer(V,x = :A,y = :B,color = :EL,Geom.rectbin),
    Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = p -> get(ColorSchemes.linear_blue_5_95_c73_n256,p),minvalue = 0,maxvalue = 200),
    Guide.xlabel("A (hours at surface)"),
    Guide.ylabel("B (timing of switch)"),
    # Guide.title("Energy for Locomotion, EL"),
    Guide.colorkey(""),
    latex_fonts
)

E = plot(
    create_contour_layers(v,:mean_E_loss,c.meta[1,:mean_E_loss])...,
    layer(examples,x = :A,y = :B,Geom.point,shape = :behavior,Theme(default_color = colorant"black",point_size = 1.5mm)),
    layer(frontier,x = :A,y = :B0,Geom.line,Theme(line_style = [:dash],default_color = colorant"gray",line_width=1mm)),
    layer(V,x = :A,y = :B,color = :mean_E_loss,Geom.rectbin),
    Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = p -> get(ColorSchemes.linear_blue_5_95_c73_n256,p)),
    Guide.xlabel("A (hours at surface)"),
    Guide.ylabel("B (timing of switch)"),
    # Guide.title("Total Energy Expenditure"),
    Guide.colorkey(""),
    latex_fonts
)

# EM = plot(
#     create_contour_layers(v,:mean_E_loss,c.meta[1,:mean_E_loss])...,
#     layer(examples,x = :A,y = :B,Geom.point,shape = :behavior,Theme(default_color = colorant"black",point_size = 1.5mm)),
#     layer(frontier,x = :A,y = :B0,Geom.line,Theme(line_style = [:dash],default_color = colorant"gray",line_width=1mm)),
#     layer(V,x = :A,y = :B,color = :mean_E_loss,Geom.rectbin),
#     Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
#     Scale.color_continuous(colormap = p -> get(ColorSchemes.linear_blue_5_95_c73_n256,p)),
#     Guide.xlabel("A (hours at surface)"),
#     Guide.ylabel("B (timing of switch)"),
#     Guide.title("Energy for Maintenance, EM"),
#     Guide.colorkey(""),
#     latex_fonts
# )

F = plot(
    create_contour_layers(v,:mean_E_gain,c.meta[1,:mean_E_gain])...,
    layer(examples,x = :A,y = :B,Geom.point,shape = :behavior,Theme(default_color = colorant"black",point_size = 1.5mm)),
    layer(frontier,x = :A,y = :B0,Geom.line,Theme(line_style = [:dash],default_color = colorant"gray",line_width=1mm)),
    layer(V,x = :A,y = :B,color = :mean_E_gain,Geom.rectbin),
    Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = p -> get(ColorSchemes.linear_blue_5_95_c73_n256,p)),
    Guide.xlabel("A (hours at surface)"),
    Guide.ylabel("B (timing of switch)"),
    # Guide.title("Mean Total Food Encountered, F"),
    Guide.colorkey(""),
    latex_fonts
)

draw(PDF(figdir*"/performance/heatmap-DE.pdf",13cm,10cm),DE)
draw(PDF(figdir*"/performance/heatmap-DE0.pdf",13cm,10cm),DE0)
draw(PDF(figdir*"/performance/heatmap-DT.pdf",13cm,10cm),DT)
draw(PDF(figdir*"/performance/heatmap-DP.pdf",13cm,10cm),DP)
draw(PDF(figdir*"/performance/heatmap-S-diurnal.pdf",13cm,10cm),S_diurnal)
draw(PDF(figdir*"/performance/heatmap-S-nearshore.pdf",13cm,10cm),S_nearshore)
draw(PDF(figdir*"/energy/heatmap-EL.pdf",13cm,10cm),EL)
draw(PDF(figdir*"/energy/heatmap-E.pdf",13cm,10cm),E)
draw(PDF(figdir*"/energy/heatmap-F.pdf",13cm,10cm),F)
