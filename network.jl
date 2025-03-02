using Random, LaTeXStrings, Plots, Statistics

# Parameters
α, β, γ, ϵ, θ, κ, ρ = 1 / 100, 0.4, 1 / 7, 1 / 3, 0.01, 1000, 0.01
N = 2000 # number of individuals
T = 400 # number of days
num_sims = 50
umax = 0.2
# Output
output_global = zeros(4, 3)
output_global_var = zeros(4, 3)
output_local = zeros(4, 3)
output_local_var = zeros(4, 3)
output_global_control = zeros(4, 4)
output_global_control_var = zeros(4, 4)
output_local_control = zeros(4, 4)
output_local_control_var = zeros(4, 4)
# Probability of being infected
function Pinf(A, status)
    return β * sum(x -> x == 2, status[A]) / length(A) #1 - (1-β)^sum(x->x==2,A)
end
# Probability of being vaccinated
function Pvax_local(A, u, status, vax)
    I = sum(x -> x == 2, status[A]) / length(A)
    v = sum(vax[A]) / length(A)
    return 1 / (1 + exp(-κ * (β * I - ρ - θ * (1 - 2 * (v + u * (1 - v))))))
end
function Pvax_global(u, status, vax)
    I = sum(x -> x == 2, status) / N
    v = sum(vax) / N
    return 1 / (1 + exp(-κ * (β * I - ρ - θ * (1 - 2 * (v + u * (1 - v))))))
end

num_edges = [5, 10, 20, 40]

for count = 1:4
    d = num_edges[count]

    ##############  No control  ##################
    control = 0 # public policy

    # Global information
    output_S = zeros(num_sims)
    output_I = zeros(num_sims)
    output_v = zeros(num_sims)
    for sim = 1:num_sims
        # Initialize variables and save initial conditions
        status = [ones(Int(N * 0.99)); 2 * ones(Int(N * 0.01))] # 1=suscpetible, 2=infected
        vax = 0.01 * ones(N)
        output_S[sim] += (sum(x -> x == 1, status) / N) / T
        output_I[sim] += (sum(x -> x == 2, status) / N) / T
        output_v[sim] += (sum(vax) / N) / T
        # Generate Erdos-Renyi random graph
        wire = d / N
        loc_neigh = [[] for i = 1:N]
        for i = 1:N-1
            for j = i+1:N
                if wire ≥ rand() # join the two players if this condition is true
                    push!(loc_neigh[i], j)
                    push!(loc_neigh[j], i)
                end
            end
        end
        for t = 1:T
            # Susceptibles become infected or vaccinated
            order = shuffle!(collect(1:N))
            for n = 1:N
                vax[n] = Pvax_global(control, status, vax)
                if !isempty(loc_neigh[order[n]])
                    # Infection
                    if status[n] == 1
                        if rand() ≤ Pinf(loc_neigh[order[n]], status)
                            status[n] = 2
                        elseif rand() ≤ 0.0001
                            status[n] = 2
                            # Vaccination
                        elseif status[n] == 1 && rand() ≤ ϵ * vax[n]
                            status[n] = 3
                        end
                    end
                end
            end

            # Recovery and resusceptibility
            for n = 1:N
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

            output_S[sim] += (sum(x -> x == 1, status) / N) / T
            output_I[sim] += (sum(x -> x == 2, status) / N) / T
            output_v[sim] += (sum(vax) / N) / T
        end
    end
    output_global[count, :] = [mean(output_S), mean(output_I), mean(output_v)]
    output_global_var[count, :] = [var(output_S, corrected=false), var(output_I, corrected=false), var(output_v, corrected=false)]

    # Local information
    output_S = zeros(num_sims)
    output_I = zeros(num_sims)
    output_v = zeros(num_sims)
    for sim = 1:num_sims
        # Initialize variables and save initial conditions
        status = [ones(Int(N * 0.99)); 2 * ones(Int(N * 0.01))] # 1=suscpetible, 2=infected
        vax = 0.01 * ones(N)
        output_S[sim] += (sum(x -> x == 1, status) / N) / T
        output_I[sim] += (sum(x -> x == 2, status) / N) / T
        output_v[sim] += (sum(vax) / N) / T
        # Generate Erdos-Renyi random graph
        wire = d / N
        loc_neigh = [[] for i = 1:N]
        for i = 1:N-1
            for j = i+1:N
                if wire ≥ rand() # join the two players if this condition is true
                    push!(loc_neigh[i], j)
                    push!(loc_neigh[j], i)
                end
            end
        end
        for t = 1:T

            # Susceptibles become infected or vaccinated
            order = shuffle!(collect(1:N))
            for n = 1:N
                if !isempty(loc_neigh[order[n]])
                    vax[n] = Pvax_local(loc_neigh[order[n]], control, status, vax)
                    # Infection
                    if status[n] == 1
                        if rand() ≤ Pinf(loc_neigh[order[n]], status)
                            status[n] = 2
                        elseif rand() ≤ 0.0001
                            status[n] = 2
                            # Vaccination
                        elseif status[n] == 1 && rand() ≤ ϵ * vax[n]
                            status[n] = 3
                        end
                    end
                end
            end

            # Recovery and resusceptibility
            for n = 1:N
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

            output_S[sim] += (sum(x -> x == 1, status) / N) / T
            output_I[sim] += (sum(x -> x == 2, status) / N) / T
            output_v[sim] += (sum(vax) / N) / T
        end
    end
    output_local[count, :] = [mean(output_S), mean(output_I), mean(output_v)]
    output_local_var[count, :] = [var(output_S, corrected=false), var(output_I, corrected=false), var(output_v, corrected=false)]

    ##############  Control  ##################

    # Global information
    output_S = zeros(num_sims)
    output_I = zeros(num_sims)
    output_v = zeros(num_sims)
    output_u = zeros(num_sims)
    for sim = 1:num_sims
        # Initialize variables and save initial conditions
        status = [ones(Int(N * 0.99)); 2 * ones(Int(N * 0.01))] # 1=suscpetible, 2=infected
        vax = 0.01 * ones(N)
        output_S[sim] += (sum(x -> x == 1, status) / N) / T
        output_I[sim] += (sum(x -> x == 2, status) / N) / T
        output_v[sim] += (sum(vax) / N) / T
        # Generate Erdos-Renyi random graph
        wire = d / N
        loc_neigh = [[] for i = 1:N]
        for i = 1:N-1
            for j = i+1:N
                if wire ≥ rand() # join the two players if this condition is true
                    push!(loc_neigh[i], j)
                    push!(loc_neigh[j], i)
                end
            end
        end
        for t = 1:T

            # Promote vaccination
            if β * sum(x -> x == 2, status) / N > ρ && sum(vax) / N < 1 / 2
                control = umax
            else
                control = 0
            end

            # Susceptibles become infected or vaccinated
            order = shuffle!(collect(1:N))
            for n = 1:N
                vax[n] = Pvax_global(control, status, vax)
                if !isempty(loc_neigh[order[n]])
                    # Infection
                    if status[n] == 1
                        if rand() ≤ Pinf(loc_neigh[order[n]], status)
                            status[n] = 2
                        elseif rand() ≤ 0.0001
                            status[n] = 2
                            # Vaccination
                        elseif status[n] == 1 && rand() ≤ ϵ * vax[n]
                            status[n] = 3
                        end
                    end
                end
            end

            # Recovery and resusceptibility
            for n = 1:N
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

            output_S[sim] += (sum(x -> x == 1, status) / N) / T
            output_I[sim] += (sum(x -> x == 2, status) / N) / T
            output_v[sim] += (sum(vax) / N) / T
            output_u[sim] += (control) / T
        end
    end
    output_global_control[count, :] = [mean(output_S), mean(output_I), mean(output_v), mean(output_u)]
    output_global_control_var[count, :] = [var(output_S, corrected=false), var(output_I, corrected=false), var(output_v, corrected=false), var(output_u, corrected=false)]

    # Local information
    output_S = zeros(num_sims)
    output_I = zeros(num_sims)
    output_v = zeros(num_sims)
    output_u = zeros(num_sims)
    for sim = 1:num_sims
        # Initialize variables and save initial conditions
        status = [ones(Int(N * 0.99)); 2 * ones(Int(N * 0.01))] # 1=suscpetible, 2=infected
        vax = 0.01 * ones(N)
        output_S[sim] += (sum(x -> x == 1, status) / N) / T
        output_I[sim] += (sum(x -> x == 2, status) / N) / T
        output_v[sim] += (sum(vax) / N) / T
        # Generate Erdos-Renyi random graph
        wire = d / N
        loc_neigh = [[] for i = 1:N]
        for i = 1:N-1
            for j = i+1:N
                if wire ≥ rand() # join the two players if this condition is true
                    push!(loc_neigh[i], j)
                    push!(loc_neigh[j], i)
                end
            end
        end
        for t = 1:T

            # Promote vaccination
            if β * sum(x -> x == 2, status) / N > ρ && sum(vax) / N < 1 / 2
                control = umax
            else
                control = 0
            end

            # Susceptibles become infected or vaccinated
            order = shuffle!(collect(1:N))
            for n = 1:N
                if !isempty(loc_neigh[order[n]])
                    vax[n] = Pvax_local(loc_neigh[order[n]], control, status, vax)
                    # Infection
                    if status[n] == 1
                        if rand() ≤ Pinf(loc_neigh[order[n]], status)
                            status[n] = 2
                        elseif rand() ≤ 0.0001
                            status[n] = 2
                            # Vaccination
                        elseif status[n] == 1 && rand() ≤ ϵ * vax[n]
                            status[n] = 3
                        end
                    end
                end
            end

            # Recovery and resusceptibility
            for n = 1:N
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

            output_S[sim] += (sum(x -> x == 1, status) / N) / T
            output_I[sim] += (sum(x -> x == 2, status) / N) / T
            output_v[sim] += (sum(vax) / N) / T
            output_u[sim] += (control) / T
        end
    end
    output_local_control[count, :] = [mean(output_S), mean(output_I), mean(output_v), mean(output_u)]
    output_local_control_var[count, :] = [var(output_S, corrected=false), var(output_I, corrected=false), var(output_v, corrected=false), var(output_u, corrected=false)]
end

##############  Plot  ##################


cpal = palette(:tab10)

plot_font = "Computer Modern"

p1 = plot(output_global[:, 1], thickness_scaling=1, linewidth=2, colour=cpal[1],
    label=L"\textrm{Global}", ylims=(-0.01, 0.51), xlims=(0.9, 4.1), grid=false, linestyle=:dash)

plot!(output_local[:, 1], thickness_scaling=1, linewidth=2, colour=cpal[2],
    label=L"\textrm{Local}", ylims=(-0.01, 0.51), xlims=(0.9, 4.1), grid=false)

plot!(output_global_control[:, 1], thickness_scaling=1, linewidth=2, colour=cpal[3],
    label=L"\textrm{Global}, \textrm{control}", ylims=(-0.01, 0.51), xlims=(0.9, 4.1), grid=false, linestyle=:dot)

plot!(output_local_control[:, 1], thickness_scaling=1, linewidth=2, colour=cpal[4],
    title=L"\textrm{Susce}\textrm{ptibles}", ylabel=L"S", xlabel=L"\delta", ylims=(-0.01, 0.51),
    label=L"\textrm{Local}, \textrm{control}", xlims=(0.9, 4.1), grid=false, legend=:outerright,
    xticks=([1, 2, 3, 4], [L"5", L"10", L"20", L"40"]), linestyle=:dashdot,
    yticks=([0, 0.25, 0.5], [L"0", L"0.25", L"0.5"]))


p2 = plot(output_global[:, 2], thickness_scaling=1, linewidth=2, colour=cpal[1],
    label=L"\textrm{Global}", ylims=(-0.001, 0.011), xlims=(0.9, 4.1), grid=false, linestyle=:dash)

plot!(output_local[:, 2], thickness_scaling=1, linewidth=2, colour=cpal[2],
    label=L"\textrm{Local}", ylims=(-0.001, 0.011), xlims=(0.9, 4.1), grid=false)

plot!(output_global_control[:, 2], thickness_scaling=1, linewidth=2, colour=cpal[3],
    label=L"\textrm{Global}, \textrm{control}", ylims=(-0.001, 0.011), xlims=(0.9, 4.1), grid=false, linestyle=:dot)

plot!(output_local_control[:, 2], thickness_scaling=1, linewidth=2, colour=cpal[4],
    title=L"\textrm{Infectious}", ylabel=L"I", xlabel=L"\delta", ylims=(-0.001, 0.011),
    label=L"\textrm{Local}, \textrm{control}", xlims=(0.9, 4.1), grid=false, legend=:outerright,
    xticks=([1, 2, 3, 4], [L"5", L"10", L"20", L"40"]), linestyle=:dashdot,
    yticks=([0, 0.005, 0.01], [L"0", L"0.05", L"0.01"]))


p3 = plot(output_global[:, 3], thickness_scaling=1, linewidth=2, colour=cpal[1],
    label=L"\textrm{Global}", ylims=(-0.01, 0.21), xlims=(0.9, 4.1), grid=false, linestyle=:dash)

plot!(output_local[:, 3], thickness_scaling=1, linewidth=2, colour=cpal[2],
    label=L"\textrm{Local}", ylims=(-0.01, 0.21), xlims=(0.9, 4.1), grid=false)

plot!(output_global_control[:, 3], thickness_scaling=1, linewidth=2, colour=cpal[3],
    label=L"\textrm{Global}, \textrm{control}", ylims=(-0.001, 0.011), xlims=(0.9, 4.1), grid=false, linestyle=:dot)

plot!(output_local_control[:, 3], thickness_scaling=1, linewidth=2, colour=cpal[4],
    title=L"\textrm{V}\textrm{accination}", ylabel=L"v", xlabel=L"\delta", ylims=(-0.01, 0.21),
    label=L"\textrm{Local}, \textrm{control}", xlims=(0.9, 4.1), grid=false, legend=:outerright,
    xticks=([1, 2, 3, 4], [L"5", L"10", L"20", L"40"]), linestyle=:dashdot,
    yticks=([0, 0.1, 0.2], [L"0", L"0.1", L"0.2"]))


p4 = plot(output_global_control[:, 4], thickness_scaling=1, linewidth=2, colour=cpal[1],
    label=L"\textrm{Global}, \textrm{control}", ylims=(-0.001, 0.021), xlims=(0.9, 4.1), grid=false, linestyle=:dash)

plot!(output_local_control[:, 4], thickness_scaling=1, linewidth=2, colour=cpal[2],
    title=L"\textrm{C}\textrm{ontrol}", ylabel=L"u", xlabel=L"\delta", ylims=(-0.001, 0.021),
    label=L"\textrm{Local}, \textrm{control}", xlims=(0.9, 4.1), grid=false, legend=:outerright,
    xticks=([1, 2, 3, 4], [L"5", L"10", L"20", L"40"]),
    yticks=([0, 0.01, 0.02], [L"0", L"0.01", L"0.02"]))

plot(p1, p2, p3, p4, layout=grid(2, 2), size=(600, 400))

savefig("Network.pdf")