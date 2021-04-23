include("requirements2.jl")
NSIM = 10000;
p = parameters2(T_LD = 24. * 25)
swims = OrderedDict(
    :Neutral => VSwim(vm_neutral2,true,[14.,1.]),
    :Binary => VSwim(vm_generic2,false,[24.,.25]),
    :DVM => VSwim(vm_generic2,false,[10.,1.]),
    :Hybrid => VSwim(vm_generic2,false,[10,.25])
);


neutral = Larvae(p,UW,pathfx2,swims[:Neutral],10^4,lite = false,no_settling = true)
binary = Larvae(p,UW,pathfx2,swims[:Binary],10^4,lite = false,no_settling = true)
dvm = Larvae(p,UW,pathfx2,swims[:DVM],10^4,lite = false,no_settling = true)
hybrid = Larvae(p,UW,pathfx2,swims[:Hybrid],10^4,lite = false,no_settling = true)

hist_df = DataFrame(
    :Neutral => neutral.sims.X[neutral.sims.t .== p[:T_LD]],
    :Binary => binary.sims.X[neutral.sims.t .== p[:T_LD]],
    :DVM => dvm.sims.X[neutral.sims.t .== p[:T_LD]],
    :Hybrid => hybrid.sims.X[neutral.sims.t .== p[:T_LD]],
    :Neutral_settle => neutral.stats.settled,
    :Binary_settle => binary.stats.settled,
    :DVM_settle => dvm.stats.settled,
    :Hybrid_settle => hybrid.stats.settled,
)

ah = plot(hist_df,x = :Neutral,Geom.histogram,color = :Neutral_settle,Scale.color_discrete_manual(my_colors...,levels = [0,1]),
    Theme(major_label_font="CMU Serif", major_label_font_size=16pt,
        minor_label_font="CMU Serif", minor_label_font_size=14pt,
        key_title_font="CMU Serif", key_title_font_size=12pt,
        key_label_font="CMU Serif", key_label_font_size=10pt,
        key_position = :none),
    # Guide.colorkey(title = "",labels = ["Delivered","Not delivered"]),
    Guide.xlabel("XT"),
    Guide.ylabel("count"),
    Coord.cartesian(xmax = 120,ymax = 800)
)
bh = plot(hist_df,x = :Binary,Geom.histogram,color = :Binary_settle,Scale.color_discrete_manual(my_colors...,levels = [0,1]),
    Theme(major_label_font="CMU Serif", major_label_font_size=16pt,
        minor_label_font="CMU Serif", minor_label_font_size=14pt,
        key_title_font="CMU Serif", key_title_font_size=12pt,
        key_label_font="CMU Serif", key_label_font_size=10pt,
        key_position = :none),
    # Guide.colorkey(title = "",labels = ["Delivered","Not delivered"]),
    Guide.xlabel("XT"),
    Guide.ylabel("count"),
    Coord.cartesian(xmax = 120,ymax = 800)
)
ch = plot(hist_df,x = :DVM,Geom.histogram,color = :DVM_settle,Scale.color_discrete_manual(my_colors...,levels = [0,1]),
    Theme(major_label_font="CMU Serif", major_label_font_size=16pt,
        minor_label_font="CMU Serif", minor_label_font_size=14pt,
        key_title_font="CMU Serif", key_title_font_size=12pt,
        key_label_font="CMU Serif", key_label_font_size=10pt,
        key_position = :none),
    # Guide.colorkey(title = "",labels = ["Delivered","Not delivered"]),
    Guide.xlabel("XT"),
    Guide.ylabel("count"),
    Coord.cartesian(xmax = 120,ymax = 800)
)
dh = plot(hist_df,x = :Hybrid,Geom.histogram,color = :Hybrid_settle,Scale.color_discrete_manual(my_colors...,levels = [0,1]),
    Theme(major_label_font="CMU Serif", major_label_font_size=16pt,
        minor_label_font="CMU Serif", minor_label_font_size=14pt,
        key_title_font="CMU Serif", key_title_font_size=12pt,
        key_label_font="CMU Serif", key_label_font_size=10pt,
        key_position = :none),
    # Guide.colorkey(title = "",labels = ["Delivered","Not delivered"]),
    Guide.xlabel("XT"),
    Guide.ylabel("count"),
    Coord.cartesian(xmax = 120,ymax = 800)
)
# draw(PDF(figdir * "/hist_neutral.pdf",13cm,10cm),ah)
draw(PDF(figdir * "/histograms/hist_binary.pdf",13cm,10cm),bh)
draw(PDF(figdir * "/histograms/hist_dvm.pdf",13cm,10cm),ch)
draw(PDF(figdir * "/histograms/hist_hybrid.pdf",13cm,10cm),dh)
