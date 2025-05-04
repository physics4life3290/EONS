#!/usr/bin/env julia
# network20.jl — 20‑species nucleosynthesis network (H → Zn)

using DifferentialEquations
using Plots

# === 1. Nucleus & Reaction types ===
struct Nucleus
    name::String
    A::Int
    Z::Int
end

struct Reaction
    reactants::Vector{Int}  # indices into species array
    products::Vector{Int}
    rate::Function          # function T9 → rate
end

# === 2. Define 20 species (index = position) ===
species = [
    Nucleus("H1",1,1), Nucleus("D",2,1), Nucleus("He3",3,2), Nucleus("He4",4,2),
    Nucleus("C12",12,6), Nucleus("N13",13,7), Nucleus("N14",14,7), Nucleus("O15",15,8), Nucleus("O16",16,8),
    Nucleus("Ne20",20,10), Nucleus("Mg24",24,12), Nucleus("Si28",28,14),
    Nucleus("S32",32,16), Nucleus("Ar36",36,18), Nucleus("Ca40",40,20),
    Nucleus("Ti44",44,22), Nucleus("Cr48",48,24), Nucleus("Fe52",52,26),
    Nucleus("Ni56",56,28), Nucleus("Zn60",60,30)
]

# === 3. Simplified rate law (dummy) ===
const T9_fixed = 3.0  # temperature in 10^9 K
dummy_rate(T9; coeff=1e-9, E=1.0) = coeff * T9^(-2/3) * exp(-E/T9^(1/3))

# === 4. Build reactions ===
reactions = Reaction[]
# H‑burning
push!(reactions, Reaction([1,1],[2], T->dummy_rate(T,coeff=1e-11,E=3.38)))   # H1+H1→D
push!(reactions, Reaction([2,1],[3], T->dummy_rate(T,coeff=1e-9, E=1.0)))    # D+H1→He3
push!(reactions, Reaction([3,3],[4,1,1], T->dummy_rate(T,coeff=1e-8, E=1.0))) # He3+He3→He4+2H

# Triple‑alpha (via C12)
push!(reactions, Reaction([4,4],[5], T->dummy_rate(T,coeff=1e-8, E=0.9)))    # He4+He4→C12

# CNO cycle
push!(reactions, Reaction([5,1],[6], T->dummy_rate(T)))  # C12+p→N13
push!(reactions, Reaction([6],[7],    T->dummy_rate(T,coeff=1e-10,E=0.5))) # N13→N14
push!(reactions, Reaction([7,1],[8], T->dummy_rate(T)))  # N14+p→O15
push!(reactions, Reaction([8],[9],    T->dummy_rate(T,coeff=1e-10,E=0.5))) # O15→O16

# α‑chain
α = 4
chain = [5,10,11,12,13,14,15,16,17,18,19,20]  # C12→…→Zn60
for (i, prod) in enumerate(chain[1:end-1])
    react = chain[i]
    nxt   = chain[i+1]
    push!(reactions, Reaction([α,react],[nxt], T->dummy_rate(T)))  # He4 + react → nxt
end

# === 5. ODE system ===
function dYdt!(dY,Y,p,t)
    T9 = p
    fill!(dY,0.0)
    for rxn in reactions
        r = rxn.rate(T9) * prod(Y[i] for i in rxn.reactants)
        for i in rxn.reactants; dY[i] -= r end
        for i in rxn.products;  dY[i] += r end
    end
end

# === 6. Initial abundances & solve ===
Y0 = zeros(Float64, length(species))
Y0[1] = 0.7   # H1
Y0[4] = 0.3   # He4

tspan = (0.0, 1e13)
prob = ODEProblem(dYdt!, Y0, tspan, T9_fixed)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-10)

# === 7. Print final abundances ===
println("Final abundances at t = $(sol.t[end]) s:")
for (i,nuc) in enumerate(species)
    println(rpad(nuc.name,6), " = ", round(sol.u[end][i], sigdigits=5))
end

# === 8. Plot evolution ===
names = [n.name for n in species]

p = plot(
    xlabel="Time (s)",
    ylabel="Abundance",
    ylims=(1e-8, 1),
    legend = :left,
    lw = 2,
    yscale = :log10
)

for i in 1:length(species)
    times = sol.t
    values = [u[i] for u in sol.u]

    # Filter out t=0.0
    ts = []
    ys = []
    for (t, y) in zip(times, values)
        if t > 0.0
            push!(ts, t)
            push!(ys, y > 0 ? y : 1e-20)  # floor zeros for log10
        end
    end

    if !isempty(ts)
        plot!(p, ts, ys, label = names[i])
    end
end

savefig(p, "network20_abundances_logy.png")
println("Plot saved to network20_abundances_logy.png")
