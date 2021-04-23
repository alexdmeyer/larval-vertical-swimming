include("requirements2.jl")
NSIM = 10000;
p = parameters2(T_LD = 24. * 5)

swims = OrderedDict(
    :Neutral => VSwim(vm_neutral2,true,[14.,1.]),
    :Binary => VSwim(vm_generic2,false,[24.,.25]),
    :DVM => VSwim(vm_generic2,false,[10.,1.]),
    :Hybrid => VSwim(vm_generic2,false,[10,.25])
);

PASSIVE = Larvae(p,UW,pathfx2,swims[:Neutral],NSIM,lite = false)
# PASSIVE = Larvae2(p,UW,pathfx2,VSwim(vm_neutral2,true,[14.,1.]),NSIM,lite = false)

ex_passive_z = plot(
    filter(trt -> (trt.n == 2),PASSIVE.sims),
    x = :t,
    y = :Z,
    color = [colorant"black"],
    Geom.line,
    # Geom.point,
    # Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    # Scale.color_continuous(colormap = p -> get(ColorSchemes.magma,p)),
    Guide.xlabel("time, t (h)"),
    Guide.ylabel("layer, Z"),
    # Guide.title("PASSIVE"),
    latex_fonts
)

ex_passive_x = plot(
    layer(
        xmin = [p[:T_PC]],
        xmax = [p[:T_LD]],
        ymin = [0],
        ymax = [p[:xh]],
        color = [my_colors[6]],
        Geom.rect
    ),
    layer(PASSIVE.sims,x = :t,y = :X,Geom.histogram2d(xbincount = p[:NT])),
    # Scale.color_log10,
    # Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Scale.color_log10(colormap = my_colormap,maxvalue = 1e3),
    Guide.xlabel("time, t (h)"),
    Guide.ylabel("offshore distance, X (km)",orientation = :vertical),
    # Guide.title("PASSIVE"),
    Coord.cartesian(ymax = 50),
    latex_fonts
)


BINARY = Larvae(p,UW,pathfx2,swims[:Binary],NSIM,lite = false)
# BINARY = Larvae2(p,UW,pathfx2,VSwim(vm_generic2,false,[24.,.25]),NSIM,lite = false)

ex_binary_z = plot(
    filter(trt -> (trt.n == 1),BINARY.sims),
    x = :t,
    y = :Z,
    color = [colorant"black"],
    Geom.line,
    # Geom.point,
    # Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    # Scale.color_continuous(colormap = p -> get(ColorSchemes.magma,p)),
    Guide.xlabel("time, t (h)"),
    Guide.ylabel("layer, Z"),
    # Guide.title("BINARY"),
    latex_fonts
)

ex_binary_x = plot(
    layer(
        xmin = [p[:T_PC]],
        xmax = [p[:T_LD]],
        ymin = [0],
        ymax = [p[:xh]],
        color = [my_colors[6]],
        Geom.rect
    ),
    layer(
        BINARY.sims,
        x = :t,
        y = :X,
        Geom.histogram2d(xbincount = p[:NT]),
    ),
    # Scale.color_log10,
    # Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Coord.cartesian(ymax = 50),
    Scale.color_log10(colormap = my_colormap,maxvalue = 1e3),
    Guide.xlabel("time, t (h)"),
    Guide.ylabel("offshore distance, X (km)",orientation = :vertical),
    # Guide.title("BINARY"),
    latex_fonts
)



DVM = Larvae(p,UW,pathfx2,swims[:DVM],NSIM,lite = false)

# DVM1 = Larvae2(p,UW,pathfx2,VSwim(vm_generic2,false,[6.,1.]),NSIM,lite = false)

ex_dvm_z = plot(
    filter(trt -> (trt.n == 1),DVM.sims),
    x = :t,
    y = :Z,
    color = [colorant"black"],
    Geom.line,
    # Geom.point,
    # Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    # Scale.color_continuous(colormap = p -> get(ColorSchemes.magma,p)),
    Guide.xlabel("time, t (h)"),
    Guide.ylabel("layer, Z"),
    # Guide.title("DVM-1"),
    latex_fonts
)

ex_dvm_x = plot(
    layer(
        xmin = [p[:T_PC]],
        xmax = [p[:T_LD]],
        ymin = [0],
        ymax = [p[:xh]],
        color = [my_colors[6]],
        Geom.rect
    ),
    layer(DVM.sims,x = :t,y = :X,Geom.histogram2d(xbincount = p[:NT])),
    # Scale.color_log10,
    # Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Coord.cartesian(ymax = 50),
    Scale.color_log10(colormap = my_colormap,maxvalue = 1e3),
    Guide.xlabel("time, t (h)"),
    Guide.ylabel("offshore distance, X (km)",orientation = :vertical),
    # Guide.title("DVM-1"),
    latex_fonts
)


HYBRID = Larvae(p,UW,pathfx2,swims[:Hybrid],NSIM,lite = false)

ex_hybrid_z = plot(
    filter(trt -> (trt.n == 1),HYBRID.sims),
    x = :t,
    y = :Z,
    color = [colorant"black"],
    Geom.line,
    # Geom.point,
    # Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    # Scale.color_continuous(colormap = p -> get(ColorSchemes.magma,p)),
    Guide.xlabel("time, t (h)"),
    Guide.ylabel("layer, Z"),
    # Guide.title("HYBRID"),
    latex_fonts
)

ex_hybrid_x = plot(
    layer(
        xmin = [p[:T_PC]],
        xmax = [p[:T_LD]],
        ymin = [0],
        ymax = [p[:xh]],
        color = [my_colors[6]],
        Geom.rect
    ),
    layer(HYBRID.sims,x = :t,y = :X,Geom.histogram2d(xbincount = p[:NT])),
    # Scale.color_log10,
    # Coord.cartesian(xmin = 0,xmax = 24,ymin = 0,ymax = 1),
    Coord.cartesian(ymax = 50),
    Scale.color_log10(colormap = my_colormap,maxvalue = 1e3),
    Guide.xlabel("time, t (h)"),
    Guide.ylabel("offshore distance, X (km)",orientation = :vertical),
    # Guide.title("HYBRID"),
    latex_fonts
)


draw(PDF(figdir * "/examples/ex-passive-z-debug.pdf",13cm,10cm),ex_passive_z)
draw(PDF(figdir * "/examples/ex-passive-x-debug.pdf",13cm,10cm),ex_passive_x)
draw(PDF(figdir * "/examples/ex-binary-z-debug.pdf",13cm,10cm),ex_binary_z)
draw(PDF(figdir * "/examples/ex-binary-x-debug.pdf",13cm,10cm),ex_binary_x)
draw(PDF(figdir * "/examples/ex-dvm-z-debug.pdf",13cm,10cm),ex_dvm_z)
draw(PDF(figdir * "/examples/ex-dvm-x-debug.pdf",13cm,10cm),ex_dvm_x)
draw(PDF(figdir * "/examples/ex-hybrid-z-debug.pdf",13cm,10cm),ex_hybrid_z)
draw(PDF(figdir * "/examples/ex-hybrid-x-debug.pdf",13cm,10cm),ex_hybrid_x)
