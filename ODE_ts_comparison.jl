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

p = [1/100,0.2,1/7,1/3,0.01,1000,0.01]
prob = ODEProblem(sipv!, u0, tspan, p)
solbeta1 = solve(prob, saveat=0.4, reltol=1e-12)

p = [1/100,0.22,1/7,1/3,0.01,1000,0.01]
prob = ODEProblem(sipv!, u0, tspan, p)
solbeta2 = solve(prob, saveat=0.4, reltol=1e-12)

p = [1/100,0.4,0.3,1/3,0.01,1000,0.01]
prob = ODEProblem(sipv!, u0, tspan, p)
solgamma1 = solve(prob, saveat=0.4, reltol=1e-12)

p = [1/100,0.4,0.32,1/3,0.01,1000,0.01]
prob = ODEProblem(sipv!, u0, tspan, p)
solgamma2 = solve(prob, saveat=0.4, reltol=1e-12)

p = [1/100,0.4,1/7,0.01,0.01,1000,0.01]
prob = ODEProblem(sipv!, u0, tspan, p)
solepsilon1 = solve(prob, saveat=0.4, reltol=1e-12)

p = [1/100,0.4,1/7,0.03,0.01,1000,0.01]
prob = ODEProblem(sipv!, u0, tspan, p)
solepsilon2 = solve(prob, saveat=0.4, reltol=1e-12)

p = [1/100,0.4,1/7,1/3,0.01,200,0.01]
prob = ODEProblem(sipv!, u0, tspan, p)
solkappa1 = solve(prob, saveat=0.4, reltol=1e-12)

p = [1/100,0.4,1/7,1/3,0.01,300,0.01]
prob = ODEProblem(sipv!, u0, tspan, p)
solkappa2 = solve(prob, saveat=0.4, reltol=1e-12)

p = [1/100,0.4,1/7,1/3,0.01,1000,0.01]
prob = ODEProblem(sipv!, u0, tspan, p)
soltheta1 = solve(prob, saveat=0.4, reltol=1e-12)

p = [1/100,0.4,1/7,1/3,0.06,1000,0.01]
prob = ODEProblem(sipv!, u0, tspan, p)
soltheta2 = solve(prob, saveat=0.4, reltol=1e-12)

p = [1/100,0.4,1/7,1/3,0.12,1000,0.01]
prob = ODEProblem(sipv!, u0, tspan, p)
soltheta3 = solve(prob, saveat=0.4, reltol=1e-12)

cpal = palette(:tab10)

plot_font = "Computer Modern"

# Beta
pbeta1 = plot(solbeta1,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash,legend=:outerright)

plot!(solbeta1,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,legend=:outerright)

plot!(solbeta1,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"\beta=0.2",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),linestyle=:dot,
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]),legend=:outerright)

pbeta2 = plot(solbeta2,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash,legend=:outerright)

plot!(solbeta2,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,legend=:outerright)

plot!(solbeta2,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"\beta=0.22",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,linestyle=:dot,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]),legend=:outerright)

# Gamma
pgamma1 = plot(solgamma1,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash,legend=:outerright)

plot!(solgamma1,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,legend=:outerright)

plot!(solgamma1,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"\gamma=0.3",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),linestyle=:dot,
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]),legend=:outerright)

pgamma2 = plot(solgamma2,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash,legend=:outerright)

plot!(solgamma2,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,legend=:outerright)

plot!(solgamma2,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"\gamma=0.32",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,linestyle=:dot,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]),legend=:outerright)

# Epsilon
pepsilon1 = plot(solepsilon1,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash,legend=:outerright)

plot!(solepsilon1,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,legend=:outerright)

plot!(solepsilon1,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"\epsilon=0.01",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),linestyle=:dot,
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]),legend=:outerright)

pepsilon2 = plot(solepsilon2,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash,legend=:outerright)

plot!(solepsilon2,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,legend=:outerright)

plot!(solepsilon2,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"\epsilon=0.03",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,linestyle=:dot,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]),legend=:outerright)

# Kappa
pkappa1 = plot(solkappa1,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash,legend=:outerright)

plot!(solkappa1,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,legend=:outerright)

plot!(solkappa1,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"\kappa=200",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),linestyle=:dot,
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]),legend=:outerright)

pkappa2 = plot(solkappa2,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash,legend=:outerright)

plot!(solkappa2,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,legend=:outerright)

plot!(solkappa2,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"\kappa=300",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,linestyle=:dot,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]),legend=:outerright)

# Theta
ptheta1 = plot(soltheta1,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash,legend=:outerright)

plot!(soltheta1,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,legend=:outerright)

plot!(soltheta1,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"\theta=0.01",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),linestyle=:dot,
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]),legend=:outerright)

ptheta2 = plot(soltheta2,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash,legend=:outerright)

plot!(soltheta2,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,legend=:outerright)

plot!(soltheta2,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"\theta=0.06",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,linestyle=:dot,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]),legend=:outerright)

ptheta3 = plot(soltheta3,idxs = [(0, 1)],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash,legend=:outerright)

plot!(soltheta3,idxs = [(0, 2)],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,legend=:outerright)

plot!(soltheta3,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"\theta=0.12",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,linestyle=:dot,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]),legend=:outerright)

plot(pbeta1,pbeta2,pgamma1,pgamma2,pepsilon1,pepsilon2,pkappa1,pkappa2, layout = grid(4, 2), size = (600,800))

savefig("param_ts.pdf")

plot(ptheta1,ptheta2,ptheta3, layout = grid(2, 2), size = (600,400))

savefig("theta_ts.pdf")