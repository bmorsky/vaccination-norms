using OptimalControl, MINPACK, Plots, NLPModelsIpopt, MadNLP

α, β, γ, ϵ, θ, κ, ρ = 1/100,0.4,1/7,1/3,0.01,1000,0.01

##### Optimal control for J1 cost functional

@def ocp1 begin
    t ∈ [0, 20], time
    x ∈ R^3, state
    S = x₁
    I = x₂
    v = x₃
    u ∈ R, control
    S(0) == 0.99
    I(0) == 0.01
    v(0) == 0.01
    ẋ(t) == [α*(1-S(t)-I(t)) - β*S(t)*I(t) - ϵ*v(t)*S(t), β*S(t)*I(t) - γ*I(t), 1/(1+exp(-κ*(β*I(t) - ρ - θ*(1-2*(v(t)+u(t)*(1-v(t))))))) - v(t)]
    0 ≤ u(t) ≤ 0.2
    0 ≤ S(t) ≤ 1
    0 ≤ I(t) ≤ 1
    0 ≤ v(t) ≤ 1
    ∫(I(t) + 0.1*u(t)) → min
end

optsol1 = solve(ocp1, acceptable_tol=1e-15)

plot(optsol1, size=(600, 450), ylims=(0,1))
savefig("optimal_control_sol1.pdf")

sol1 = zeros(400,3)
for n=1:400
    sol1[n,1] = state(optsol1)(20*n/400)[2]
    sol1[n,2] = state(optsol1)(20*n/400)[3]
    sol1[n,3] = control(optsol1)(20*n/400)
end

##### Optimal control for J2 cost functional

@def ocp2 begin
    t ∈ [0, 20], time
    x ∈ R^3, state
    S = x₁
    I = x₂
    v = x₃
    u ∈ R, control
    S(0) == 0.99
    I(0) == 0.01
    v(0) == 0.01
    ẋ(t) == [α*(1-S(t)-I(t)) - β*S(t)*I(t) - ϵ*v(t)*S(t), β*S(t)*I(t) - γ*I(t), 1/(1+exp(-κ*(β*I(t) - ρ - θ*(1-2*(v(t)+u(t)*(1-v(t))))))) - v(t)]
    0 ≤ u(t) ≤ 0.2
    0 ≤ S(t) ≤ 1
    0 ≤ I(t) ≤ 1
    0 ≤ v(t) ≤ 1
    ∫(40*I(t)^2 + 0.1*u(t)) → min
end

optsol2 = solve(ocp2, acceptable_tol=1e-15)

plot(optsol2, size=(600, 450), ylims=(0,1))
savefig("optimal_control_sol2.pdf")

sol3 = zeros(400,3)
for n=1:400
    sol3[n,1] = state(optsol2)(20*n/400)[2]
    sol3[n,2] = state(optsol2)(20*n/400)[3]
    sol3[n,3] = control(optsol2)(20*n/400)
end

using DifferentialEquations, LaTeXStrings

##### Non-control case
function sipv!(du, u, p, t)
    α, β, γ, ϵ, θ, κ, ρ = p
    S, I, v = u
    du[1] = α*(1-S-I) - β*S*I - ϵ*v*S
    du[2] = β*S*I - γ*I
    du[3] = 1/(1+exp(-κ*(β*I - ρ - θ*(1-2*v)))) - v
end

u0 = [0.99; 0.01; 0.01]
tspan = (0.0, 20.0)
p = [1/100,0.4,1/7,1/3,0.01,1000,0.01]

prob = ODEProblem(sipv!, u0, tspan, p)
sol2 = solve(prob, saveat=0.1, reltol=1e-12)

cpal = palette(:tab10)

plot_font = "Computer Modern"

xvals = 20*collect(1:400)/400
ruleofthumb = ρ*ones(400)/β

p1 = plot(xvals,sol1[:,1],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"I, \textrm{with} \: \textrm{control}",ylims=(-0.01,0.21),xlims=(-0.1,10),grid=false,linestyle=:dash)

plot!(xvals,sol1[:,2],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"v, \textrm{with} \: \textrm{control}",ylims=(-0.01,0.21),xlims=(-0.1,10),grid=false)

plot!(xvals,sol1[:,3],thickness_scaling = 1,linewidth=2,colour=cpal[3],
label=L"u, \textrm{with} \: \textrm{control}",ylims=(-0.01,0.21),xlims=(-0.1,10),grid=false,linestyle=:dot)

plot!(sol2,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[4],
title=L"\mathcal{J}_1",ylabel=L"I,v,u",xlabel=L"t",ylims=(-0.01,0.21),
label=L"v, \textrm{without} \: \textrm{control}",xlims=(-0.1,10),grid=false,legend=:outerright,
xticks = ([0,2.5,5,7.5,10],[L"0",L"2.5", L"5", L"7.5",L"10"]),linestyle=:dashdot,
yticks = ([0,ρ/β,0.1,0.2],[L"0", L"ρ/β", L"0.1", L"0.2"]))

plot!(xvals,ruleofthumb,thickness_scaling = 1,linewidth=2,colour=cpal[5],
label=L"ρ/β",ylims=(-0.01,0.21),xlims=(-0.1,10),grid=false)

p2 = plot(xvals,sol3[:,1],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"I, \textrm{with} \: \textrm{control}",ylims=(-0.01,0.21),xlims=(-0.1,10),grid=false,linestyle=:dash)

plot!(xvals,sol3[:,2],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"v, \textrm{with} \: \textrm{control}",ylims=(-0.01,0.21),xlims=(-0.1,10),grid=false)

plot!(xvals,sol3[:,3],thickness_scaling = 1,linewidth=2,colour=cpal[3],
label=L"u, \textrm{with} \: \textrm{control}",ylims=(-0.01,0.21),xlims=(-0.1,10),grid=false,linestyle=:dot)

plot!(sol2,idxs = [(0, 3)],thickness_scaling = 1,linewidth=2,colour=cpal[4],
title=L"\mathcal{J}_2",ylabel=L"I,v,u",xlabel=L"t",ylims=(-0.01,0.21),
label=L"v, \textrm{without} \: \textrm{control}",xlims=(-0.1,10),grid=false,legend=:outerright,
xticks = ([0,2.5,5,7.5,10],[L"0",L"2.5", L"5", L"7.5",L"10"]),linestyle=:dashdot,
yticks = ([0,ρ/β,0.1,0.2],[L"0", L"ρ/β", L"0.1", L"0.2"]))

plot!(xvals,ruleofthumb,thickness_scaling = 1,linewidth=2,colour=cpal[5],
label=L"ρ/β",ylims=(-0.01,0.21),xlims=(-0.1,10),grid=false)

plot(p1,p2, layout = grid(2,1), size = (600,400))

savefig("optimal_control.pdf")