using Random,LaTeXStrings,Plots

# Parameters
α, β, γ, ϵ, θ, κ, ρ = 1/100,0.4,1/7,1/3,0.01,1000,0.01
N = 2000 # number of individuals
T = 400 # number of days
umax = 0.2
# Output
output_local_d5 = zeros(T+1,3)
output_local_d40 = zeros(T+1,3)
output_local_d5_control = zeros(T+1,4)
output_local_d40_control = zeros(T+1,4)
# Probability of being infected
function Pinf(A)
    return β*sum(x->x==2,status[A])/length(A) #1 - (1-β)^sum(x->x==2,A)
end
# Probability of being vaccinated
function Pvax_local(A,u)
    I = sum(x->x==2,status[A])/length(A)
    v = sum(vax[A])/length(A)
    return 1/(1+exp(-κ*(β*I - ρ - θ*(1-2*(v+u*(1-v))))))
end

##############  d=5  ##################
control = 0 # public policy
# Generate Erdos-Renyi random graph
d = 5
wire = d/N
loc_neigh = [[] for i = 1:N]
for i=1:N-1
    for j=i+1:N
        if wire ≥ rand() # join the two players if this condition is true
            push!(loc_neigh[i],j)
            push!(loc_neigh[j],i)
        end
    end
end
# Local information no control
# Initialize variables and save initial conditions
status = [ones(Int(N*0.99)); 2*ones(Int(N*0.01))] # 1=suscpetible, 2=infected
vax = 0.01*ones(N)
output_local_d5[1,:] = [sum(x->x==1,status)/N, sum(x->x==2,status)/N, sum(vax)/N]
for t=1:T

    # Susceptibles become infected or vaccinated
    order = shuffle!(collect(1:N))
    for n=1:N
        if !isempty(loc_neigh[order[n]])
            vax[n] = Pvax_local(loc_neigh[order[n]],control)
            # Infection
            if status[n] == 1
                if rand() ≤ Pinf(loc_neigh[order[n]])
                    status[n] = 2
                elseif rand() ≤ 0.0001
                    status[n] = 2
                # Vaccination
                elseif status[n] == 1 && rand() ≤ ϵ*vax[n]
                    status[n] = 3
                end
            end
        end
    end

    # Recovery and resusceptibility
    for n=1:N
        if status[n] != 1
            if status[n] == 2
                if rand() ≤ γ
                    status[n] = 3
                end
            elseif status[n] == 3
                if rand() ≤ α
                    status[n] = 1
                end
            end
        end
    end

    output_local_d5[t+1,:] = [sum(x->x==1,status)/N, sum(x->x==2,status)/N, sum(vax)/N]
end
# Local information control
# Initialize variables and save initial conditions
status = [ones(Int(N*0.99)); 2*ones(Int(N*0.01))] # 1=suscpetible, 2=infected
vax = 0.01*ones(N)
output_local_d5_control[1,:] = [sum(x->x==1,status)/N, sum(x->x==2,status)/N, sum(vax)/N, 0]
for t=1:T

    # Promote vaccination
    if β*sum(x->x==2,status)/N > ρ && sum(vax)/N < 1/2
        control = umax
    else
        control = 0
    end

    # Susceptibles become infected or vaccinated
    order = shuffle!(collect(1:N))
    for n=1:N
        if !isempty(loc_neigh[order[n]])
            vax[n] = Pvax_local(loc_neigh[order[n]],control)
            # Infection
            if status[n] == 1
                if rand() ≤ Pinf(loc_neigh[order[n]])
                    status[n] = 2
                elseif rand() ≤ 0.0001
                    status[n] = 2
                # Vaccination
                elseif status[n] == 1 && rand() ≤ ϵ*vax[n]
                    status[n] = 3
                end
            end
        end
    end

    # Recovery and resusceptibility
    for n=1:N
        if status[n] != 1
            if status[n] == 2
                if rand() ≤ γ
                    status[n] = 3
                end
            elseif status[n] == 3
                if rand() ≤ α
                    status[n] = 1
                end
            end
        end
    end

    output_local_d5_control[t+1,:] = [sum(x->x==1,status)/N, sum(x->x==2,status)/N, sum(vax)/N, control]
end

##############  d=40  ##################
control = 0 # public policy
# Generate Erdos-Renyi random graph
d = 40
wire = d/N
loc_neigh = [[] for i = 1:N]
for i=1:N-1
    for j=i+1:N
        if wire ≥ rand() # join the two players if this condition is true
            push!(loc_neigh[i],j)
            push!(loc_neigh[j],i)
        end
    end
end
# Local information no control
# Initialize variables and save initial conditions
status = [ones(Int(N*0.99)); 2*ones(Int(N*0.01))] # 1=suscpetible, 2=infected
vax = 0.01*ones(N)
output_local_d40[1,:] = [sum(x->x==1,status)/N, sum(x->x==2,status)/N, sum(vax)/N]
for t=1:T

    # Susceptibles become infected or vaccinated
    order = shuffle!(collect(1:N))
    for n=1:N
        if !isempty(loc_neigh[order[n]])
            vax[n] = Pvax_local(loc_neigh[order[n]],control)
            # Infection
            if status[n] == 1
                if rand() ≤ Pinf(loc_neigh[order[n]])
                    status[n] = 2
                elseif rand() ≤ 0.0001
                    status[n] = 2
                # Vaccination
                elseif status[n] == 1 && rand() ≤ ϵ*vax[n]
                    status[n] = 3
                end
            end
        end
    end

    # Recovery and resusceptibility
    for n=1:N
        if status[n] != 1
            if status[n] == 2
                if rand() ≤ γ
                    status[n] = 3
                end
            elseif status[n] == 3
                if rand() ≤ α
                    status[n] = 1
                end
            end
        end
    end

    output_local_d40[t+1,:] = [sum(x->x==1,status)/N, sum(x->x==2,status)/N, sum(vax)/N]
end
# Local information
# Initialize variables and save initial conditions
status = [ones(Int(N*0.99)); 2*ones(Int(N*0.01))] # 1=suscpetible, 2=infected
vax = 0.01*ones(N)
output_local_d40_control[1,:] = [sum(x->x==1,status)/N, sum(x->x==2,status)/N, sum(vax)/N, 0]
for t=1:T

    # Promote vaccination
    if β*sum(x->x==2,status)/N > ρ && sum(vax)/N < 1/2
        control = umax
    else
        control = 0
    end

    # Susceptibles become infected or vaccinated
    order = shuffle!(collect(1:N))
    for n=1:N
        if !isempty(loc_neigh[order[n]])
            vax[n] = Pvax_local(loc_neigh[order[n]],control)
            # Infection
            if status[n] == 1
                if rand() ≤ Pinf(loc_neigh[order[n]])
                    status[n] = 2
                elseif rand() ≤ 0.0001
                    status[n] = 2
                # Vaccination
                elseif status[n] == 1 && rand() ≤ ϵ*vax[n]
                    status[n] = 3
                end
            end
        end
    end

    # Recovery and resusceptibility
    for n=1:N
        if status[n] != 1
            if status[n] == 2
                if rand() ≤ γ
                    status[n] = 3
                end
            elseif status[n] == 3
                if rand() ≤ α
                    status[n] = 1
                end
            end
        end
    end

    output_local_d40_control[t+1,:] = [sum(x->x==1,status)/N, sum(x->x==2,status)/N, sum(vax)/N, control]
end


##############  Plot  ##################


cpal = palette(:tab10)

plot_font = "Computer Modern"

# Local, no control, d=5
p1=plot(output_local_d5[:,1],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash)

plot!(output_local_d5[:,2],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false)

plot!(output_local_d5[:,3],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"\delta=5",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,legend=:outerright,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),linestyle=:dot,
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]))

# Local, no control, d=40
p2=plot(output_local_d40[:,1],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash)

plot!(output_local_d40[:,2],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false)

plot!(output_local_d40[:,3],thickness_scaling = 1,linewidth=2,colour=cpal[3],
title=L"\delta=40",ylabel=L"S,I,v",xlabel=L"t",ylims=(-0.01,1.01),
label=L"v",xlims=(-1,400),grid=false,legend=:outerright,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),linestyle=:dashdot,
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]))

# Local, control, d=5
p3=plot(output_local_d5_control[:,1],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash)

plot!(output_local_d5_control[:,2],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false)

plot!(output_local_d5_control[:,3],thickness_scaling = 1,linewidth=2,colour=cpal[3],
label=L"v",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dot)

plot!(output_local_d5_control[:,4],thickness_scaling = 1,linewidth=2,colour=cpal[4],
title=L"\delta=5, \textrm{control}",ylabel=L"S,I,v,u",xlabel=L"t",ylims=(-0.01,1.01),
label=L"u",xlims=(-1,400),grid=false,legend=:outerright,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),linestyle=:dashdot,
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]))

# Local, control, d=40
p4=plot(output_local_d40_control[:,1],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"S",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dash)

plot!(output_local_d40_control[:,2],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"I",ylims=(-0.01,1.01),xlims=(-1,400),grid=false)

plot!(output_local_d40_control[:,3],thickness_scaling = 1,linewidth=2,colour=cpal[3],
label=L"v",ylims=(-0.01,1.01),xlims=(-1,400),grid=false,linestyle=:dot)

plot!(output_local_d40_control[:,4],thickness_scaling = 1,linewidth=2,colour=cpal[4],
title=L"\delta=40, \textrm{control}",ylabel=L"S,I,v,u",xlabel=L"t",ylims=(-0.01,1.01),
label=L"u",xlims=(-1,400),grid=false,legend=:outerright,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),linestyle=:dashdot,
yticks = ([0,0.5,1],[L"0", L"0.5", L"1"]))

# Infections
p5=plot(output_local_d5[:,2],thickness_scaling = 1,linewidth=2,colour=cpal[1],
label=L"\delta=5",ylims=(-0.001,0.081),xlims=(-1,400),grid=false,linestyle=:dash)

plot!(output_local_d40[:,2],thickness_scaling = 1,linewidth=2,colour=cpal[2],
label=L"\delta=40",ylims=(-0.001,0.081),xlims=(-1,400),grid=false)

plot!(output_local_d5_control[:,2],thickness_scaling = 1,linewidth=2,colour=cpal[3],
label=L"\delta=5, \textrm{control}",ylims=(-0.001,0.081),xlims=(-1,400),grid=false,linestyle=:dot)

plot!(output_local_d40_control[:,2],thickness_scaling = 1,linewidth=2,colour=cpal[4],
title=L"\textrm{Infectious}",ylabel=L"I",xlabel=L"t",ylims=(-0.001,0.081),
label=L"\delta=40, \textrm{control}",xlims=(-1,400),grid=false,legend=:outerright,
xticks = ([0,100,200,300,400],[L"0",L"100", L"200", L"300",L"400"]),linestyle=:dashdot,
yticks = ([0,0.04,0.08],[L"0", L"0.04", L"0.08"]))

nullplot = plot(legend=false,grid=false,foreground_color=:white)

l = @layout([
    p1  p2;
    p3  p4;
    p5
])

plot(p1,p2,p3,p4,p5, layout = l, size = (600,600))

savefig("Network_local_ts.pdf")
