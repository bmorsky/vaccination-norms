using DifferentialEquations, LaTeXStrings, Plots

T=5000.0

outalpha = zeros(100,2)
outbeta = zeros(100,2)
outgamma = zeros(100,2)
outepsilon = zeros(100,2)
outtheta = zeros(100,2)
outkappa = zeros(100,2)
outrho = zeros(100,2)

function sipv!(du, u, p, t)
    α, β, γ, ϵ, θ, κ, rho = p
    S, I, v = u
    du[1] = α*(1-S-I) - β*S*I - ϵ*v*S
    du[2] = β*S*I - γ*I
    du[3] = 1/(1+exp(-κ*(β*I - rho - θ*(1-2*v)))) - v
end

function avgtraj(p)
    u0 = [0.99; 0.01; 0.01]
    tspan = (0.0, T)
    prob = ODEProblem(sipv!, u0, tspan, p)
    sol = solve(prob,saveat=0.1, reltol=1e-12, abstol=1e-12)
    return sum(sol[2,:])/(T*10)
end

for i=1:100
    # default p = [α,β,γ,ϵ,θ,κ,ρ] = [1/100,0.4,1/7,1/3,0.01,1000,0.01]
    param = i/100
    outalpha[i,:] = [param/20, avgtraj([param/20,0.4,1/7,1/3,0.01,1000,0.01])]
    outbeta[i,:] = [param, avgtraj([1/100,param,1/7,1/3,0.01,1000,0.01])]
    outgamma[i,:] = [param, avgtraj([1/100,0.4,param,1/3,0.01,1000,0.01])]
    outepsilon[i,:] = [param, avgtraj([1/100,0.4,1/7,param,0.01,1000,0.01])]
    outtheta[i,:] = [param*0.12, avgtraj([1/100,0.4,1/7,1/3,param*0.12,1000,0.01])]
    outkappa[i,:] = [1000*param, avgtraj([1/100,0.4,1/7,1/3,0.01,1000*param,0.01])]
    outrho[i,:] = [param/20, avgtraj([1/100,0.4,1/7,1/3,0.01,1000,param/20])]
end

cpal = palette(:tab10)

p1=plot(outalpha[:,1],outalpha[:,2],
thickness_scaling = 1,linewidth=2,legend=false,
ylabel=L"\bar{I}(5000)", ylims=(-0.001,0.025), yticks = ([0,0.01,0.02],[L"0", L"0.01", L"0.02"]),
xlims=(-0.004,0.04), xticks = ([0,0.02,0.04],[L"0", L"0.02", L"0.04"]),
grid=false,colour=cpal[1],
xlabel=L"\alpha")

p2=plot(outbeta[:,1],outbeta[:,2],
thickness_scaling = 1,linewidth=2,legend=false,
ylabel=L"\bar{I}(5000)", ylims=(-0.005,0.025), yticks = ([0,0.01,0.02],[L"0", L"0.01", L"0.02"]),
xlims=(-0.1,1), xticks = ([0,0.5,1],[L"0", L"0.5", L"1"]),
grid=false,colour=cpal[2],
xlabel=L"\beta")

p3=plot(outgamma[:,1],outgamma[:,2],
thickness_scaling = 1,linewidth=2,legend=false,
ylabel=L"\bar{I}(5000)", ylims=(-0.01,0.061), yticks = ([0,0.03,0.06],[L"0", L"0.03", L"0.06"]),
xlims=(-0.1,1), xticks = ([0,0.5,1],[L"0", L"0.5", L"1"]),
grid=false,colour=cpal[3],
xlabel=L"\gamma")

p4=plot(outepsilon[:,1],outepsilon[:,2],
thickness_scaling = 1,linewidth=2,legend=false,
ylabel=L"\bar{I}(5000)", ylims=(-0.005,0.025), yticks = ([0,0.01,0.02],[L"0", L"0.01", L"0.02"]),
xlims=(-0.1,1), xticks = ([0,0.5,1],[L"0", L"0.5", L"1"]),
grid=false,colour=cpal[4],
xlabel=L"\epsilon")

p5=plot(outtheta[:,1],outtheta[:,2],
thickness_scaling = 1,linewidth=2,legend=false,
ylabel=L"\bar{I}(5000)", ylims=(-0.01,0.061), yticks = ([0,0.03,0.06],[L"0", L"0.03", L"0.06"]),
xlims=(-0.012,0.12), xticks = ([0,0.06,0.12],[L"0", L"0.06", L"0.12"]),
grid=false,colour=cpal[5],
xlabel=L"\theta")

p6=plot(outkappa[:,1],outkappa[:,2],
thickness_scaling = 1,linewidth=2,legend=false,
ylabel=L"\bar{I}(5000)", ylims=(-0.005,0.025), yticks = ([0,0.01,0.02],[L"0", L"0.01", L"0.02"]),
xlims=(-100,1000),xticks = ([0,500,1000],[L"0", L"500", L"1000"]),
grid=false,colour=cpal[6],
xlabel=L"\kappa")

p7=plot(outrho[:,1],outrho[:,2],
thickness_scaling = 1,linewidth=2,legend=false,
ylabel=L"\bar{I}(5000)", ylims=(-0.01,0.06), yticks = ([0,0.03,0.06],[L"0", L"0.03", L"0.06"]),
xlims=(-0.004,0.04),xticks = ([0,0.02,0.04],[L"0", L"0.03", L"0.06"]),
grid=false,colour=cpal[7],
xlabel=L"\rho")

plot(p1, p2, p3, p4, p5, p6, p7, layout = grid(3, 3), size = (500,300))

savefig("param.pdf")