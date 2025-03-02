using DifferentialEquations, LaTeXStrings, Plots

function sipv!(du, u, p, t)
    α, β, γ, ϵ, θ, κ, ρ = p
    S, I, v = u
    du[1] = α*(1-S-I) - β*S*I - ϵ*v*S
    du[2] = β*S*I - γ*I
    du[3] = 1/(1+exp(-κ*(β*I - ρ - θ*(1-2*v)))) - v
end

tspan = (0.0, 400.0)

p = [1/100,0.4,1/7,1/3,0.12,1000,0.01]
u0 = [0.99; 0.01; 0.01]
prob = ODEProblem(sipv!, u0, tspan, p)
sol1 = solve(prob, saveat=0.4, reltol=1e-12)

p = [1/100,0.4,1/7,1/3,0.12,1000,0.01]
u0 = [0.99; 0.01; 0.99]
prob = ODEProblem(sipv!, u0, tspan, p)
sol2 = solve(prob, saveat=0.4, reltol=1e-12)

cpal = palette(:tab10)

plot_font = "Computer Modern"

p1 = plot(sol1,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash)

plot!(sol1,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false)

plot!(sol1,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"v_0=0.01",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,legend=:outerright,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),linestyle=:dot,
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]))

p2 = plot(sol2,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash)

plot!(sol2,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false)

plot!(sol2,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"v_0=0.99",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,legend=:outerright,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),linestyle=:dot,
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]))

nullplot = plot(legend=false,grid=false,foreground_color=:white)

plot(p1,p2, grid=(1,2), size = (600,200))

savefig("endemic.pdf")