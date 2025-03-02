using DifferentialEquations, LaTeXStrings, Plots

function sipv!(du, u, p, t)
    α, β, γ, ϵ, θ, κ, ρ = p
    S, I, v = u
    du[1] = α*(1-S-I) - β*S*I - ϵ*v*S
    du[2] = β*S*I - γ*I
    du[3] = 1/(1+exp(-κ*(β*I - ρ - θ*(1-2*v)))) - v
end

u0 = [0.99; 0.01; 0.01]
tspan = (0.0, 400.0)

p = [1/100,0.4,1/7,1/3,0.01,1000,0.01]
prob = ODEProblem(sipv!, u0, tspan, p)
sol1 = solve(prob, saveat=0.4, reltol=1e-12)

p = [1/100,0.4,1/7,0,0.01,1000,0.01]
prob = ODEProblem(sipv!, u0, tspan, p)
sol2 = solve(prob, saveat=0.4, reltol=1e-12)

p = [1/100,0.4,1/7,1/3,0.0,1000,0.01]
prob = ODEProblem(sipv!, u0, tspan, p)
sol3 = solve(prob, saveat=0.4, reltol=1e-12)

cpal = palette(:tab10)

plot_font = "Computer Modern"

p1 = plot(sol1,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash)

plot!(sol1,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false)

plot!(sol1,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"1/\epsilon = 3, \theta=0.01",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,legend=:outerright,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),linestyle=:dot,
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]))

p2 = plot(sol2,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash)

plot!(sol2,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
title=L"\textrm{N}\textrm{o} \: \textrm{vaccine}",ylabel=L"S,I",xlabel=L"t",ylims=(-0.01,1.01),
label=L"I",xlims=(-1,400),grid=false,legend=:outerright,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]))

p3 = plot(sol3,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash)

plot!(sol3,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false)

plot!(sol3,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"\textrm{N}\textrm{o} \: \textrm{social}\: \textrm{norm}",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,legend=:outerright,linestyle=:dot,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]))

p4= plot(sol1,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"1/\epsilon = 3, \theta=0.01",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash,right_margin=50*Plots.mm)

plot!(sol2,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"\textrm{N}\textrm{o} \: \textrm{vaccine}",ylims=(-0.01,1.01),xlims=(-1,400),grid=false)

plot!(sol3,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title="",ylabel=L"I",xlabel=L"t",ylims=(-0.004,0.4004),
label=L"\textrm{N}\textrm{o} \: \textrm{social}\: \textrm{norm}",
xlims=(-1,400),grid=false,legend=:topright,linestyle=:dot,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),
yticks = ([0,0.2,0.4],[L"0", L"0.2", L"0.4"]))

nullplot = plot(legend=false,grid=false,foreground_color=:white)

l = @layout([
    p1  p2;
    p3{0.885w} p4{0.735w} nullplot{0.255w} # e is an invisible white space column
])

plot(p1,p2,p3,p4, layout = l, size = (600,400))

savefig("cycles.pdf")