### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 4c9c24cc-b865-4825-a841-f717120d27d2
begin
		using Pkg
		Pkg.activate(".")
		using Colors
		using AlgebraicPetri
		using Catlab
		using GraphViz
		using JSON
		using DataFrames
		using CSV
		using XLSX
end

# ╔═╡ 32c8703f-6aa3-46be-a91b-ff36225d6bd8
module EnzymeReactions

using AlgebraicPetri
using Catlab.Programs
using Catlab.Graphics
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Distributions

using DifferentialEquations
using Plots
using Plots.PlotMeasures
export ob, ode,
       inactivate, bindunbind, degrade,
       enzX, enzXY, enzXsubY,
       enz, enz_enz, enz_sub,
       enzyme_uwd, enzyme_generators

sep_sym = "↦"
ob(type, x) = codom(Open([first(x)], LabelledReactionNet{type,Number}(x), [first(x)])).ob;
ob(x) = codom(Open([x], LabelledPetriNet(x), [x])).ob;

ode(x::Union{AbstractReactionNet{Distribution, Number},AbstractLabelledReactionNet{Distribution, Number}}, t) = begin
  β = mean.(rates(x))
  ODEProblem(vectorfield(x), concentrations(x), t, β)
end
ode(x, t) = ODEProblem(vectorfield(x), concentrations(x), t, rates(x));

function inactivate(in,on::T) where T
  inact = Symbol(first(in), :_inact)
  Open(LabelledReactionNet{T,Number}(unique((in, inact=>0)), ((Symbol(:inact_,first(in)),on),first(in)=>inact)))
end;

function bindunbind(in1, in2, on::T, off::T) where T
  out = Symbol(first(in1), sep_sym,first(in2))
  Open(LabelledReactionNet{T,Number}(unique((in1, in2,out=>0)), ((Symbol(:bind_,first(in1), sep_sym,first(in2)),on),(first(in1),first(in2))=>out),
                                                                ((Symbol(:unbind_,out),off),out=>(first(in1),first(in2)))))
end;

function degrade(prod1,prod2,on::T) where T
  in = Symbol(first(prod1), sep_sym,first(prod2))
  prod2str = String(first(prod2))
  degprod2 = Symbol(endswith(prod2str, "inact") ? first(prod2str) : prod2str, :_deg)
  Open(LabelledReactionNet{T,Number}(unique((in=>0, prod1,degprod2=>0)), ((Symbol(:deg_,in),on),in=>(first(prod1),degprod2))));
end;

function inactivate(in)
  inact = Symbol(in, :_inact)
  Open(LabelledPetriNet(unique((in, inact)), (Symbol(:inact_,in),in=>inact)))
end;

function bindunbind(in1, in2)
  out = Symbol(in1, sep_sym,in2)
  Open(LabelledPetriNet(unique((in1, in2,out)), (Symbol(:bind_,in1, sep_sym,in2),(in1,in2)=>out),
                                                (Symbol(:unbind_,out),out=>(in1,in2))))
end;

function degrade(prod1,prod2)
  in = Symbol(prod1, sep_sym,prod2)
  prod2str = String(prod2)
  degprod2 = Symbol(endswith(prod2str, "inact") ? first(prod2str) : prod2str, :_deg)
  Open(LabelledPetriNet(unique((in, prod1,degprod2)), (Symbol(:deg_,in),in=>(prod1,degprod2))));
end;

# ## Cathepsin *X* reacting with itself

enzX = @relation (X, Xinact, Xdeg) where (X, Xinact, Xdeg, XX, XXinact) begin
  inactX(X, Xinact)
  bindXX(X, XX)
  degXX(XX, X, Xdeg)
  bindXXinact(X, Xinact, XXinact)
  degXXinact(XXinact, X, Xdeg)
end

# ## Cathepsin *X* reacting with Substrate *Y*

enzXsubY = @relation (X, Xinact, Xdeg, Y, Ydeg) where (X, Xinact, Xdeg, Y, XY, Ydeg) begin
  bindXY(X, Y, XY)
  degXY(XY, X, Ydeg)
end

# ## Cathepsin *X* reacting with Cathepsin *Y*

enzXY = @relation (X, Xinact, Xdeg, Y, Yinact, Ydeg) where (X, Xinact, Xdeg, Y, Yinact, Ydeg, XY, XYinact) begin
  bindXY(X, Y, XY)
  degXY(XY, X, Ydeg)
  bindXYinact(X, Yinact, XYinact)
  degXYinact(XYinact, X, Ydeg)
end

function enz(rxns, cat)
  catsym = first(cat)
  obtype = valtype(rates(apex(first(last(first(rxns))))))
  out = oapply(enzX, Dict([:inactX, :bindXX, :degXX, :bindXXinact, :degXXinact] .=> rxns[catsym]), Dict(
    :X=>ob(obtype, cat),
    :Xinact=>ob(obtype, Symbol(catsym,:_inact)=>0),
    :Xdeg=>ob(obtype, Symbol(catsym,:_deg)=>0),
    :XX=>ob(obtype, Symbol(catsym, sep_sym,catsym)=>0),
    :XXinact=>ob(obtype, Symbol(catsym, sep_sym,catsym,:_inact)=>0)))
  bundle_legs(out, [[1,2,3]])
end

function enz_sub(rxns, cat1, sub)
  catsym = first(cat1)
  subsym = first(sub)
  catsub = Symbol(catsym, sep_sym, subsym)
  obtype = valtype(rates(apex(first(last(first(rxns))))))
  out = oapply(enzXsubY, Dict([:bindXY, :degXY] .=> rxns[catsub]), Dict(
    :X=>ob(obtype, cat1),
    :Xinact=>ob(obtype, Symbol(catsym,:_inact)=>0),
    :Xdeg=>ob(obtype, Symbol(catsym,:_deg)=>0),
    :Y=>ob(obtype, sub),
    :XY=>ob(obtype, Symbol(catsym,sep_sym,subsym)=>0),
    :Ydeg=>ob(obtype, Symbol(subsym,:_deg)=>0)))
  bundle_legs(out, [[1,2,3], [4,5]])
end

function enz_enz(rxns, cat1, cat2)
  cat1sym = first(cat1)
  cat2sym = first(cat2)
  catcat = Symbol(cat1sym, sep_sym, cat2sym)
  obtype = valtype(rates(apex(first(last(first(rxns))))))
  out = oapply(enzXY, Dict([:bindXY, :degXY, :bindXYinact, :degXYinact] .=> rxns[catcat]), Dict(
    :X=>ob(obtype, cat1),
    :Xinact=>ob(obtype, Symbol(cat1sym,:_inact)=>0),
    :Xdeg=>ob(obtype, Symbol(cat1sym,:_deg)=>0),
    :Y=>ob(obtype, cat2),
    :Yinact=>ob(obtype, Symbol(cat2sym,:_inact)=>0),
    :Ydeg=>ob(obtype, Symbol(cat2sym,:_deg)=>0),
    :XY=>ob(obtype, catcat=>0),
    :XYinact=>ob(obtype, Symbol(catcat,:_inact)=>0)))
  bundle_legs(out, [[1,2,3], [4,5,6]])
end

function enz(cat::Symbol)
  catsym = cat
  out = oapply(enzX, Dict(:inactX=>inactivate(cat), :bindXX=>bindunbind(cat, cat), :degXX=>degrade(cat, cat),
                          :bindXXinact=>bindunbind(cat, Symbol(cat,:_inact)),
                          :degXXinact=>degrade(cat, Symbol(cat, :_inact))), Dict(
    :X=>ob(cat),
    :Xinact=>ob(Symbol(catsym,:_inact)),
    :Xdeg=>ob(Symbol(catsym,:_deg)),
    :XX=>ob(Symbol(catsym, sep_sym,catsym)),
    :XXinact=>ob(Symbol(catsym, sep_sym,catsym,:_inact))))
  bundle_legs(out, [[1,2,3]])
end

function enz_sub(cat1::Symbol, sub::Symbol)
  catsym = cat1
  subsym = sub
  catsub = Symbol(catsym, sep_sym, subsym)
  out = oapply(enzXsubY, Dict(:bindXY=>bindunbind(cat1, sub), :degXY=>degrade(cat1, sub)), Dict(
    :X=>ob(cat1),
    :Xinact=>ob(Symbol(catsym,:_inact)),
    :Xdeg=>ob(Symbol(catsym,:_deg)),
    :Y=>ob(sub),
    :XY=>ob(Symbol(catsym, sep_sym,subsym)),
    :Ydeg=>ob(Symbol(subsym,:_deg))))
  bundle_legs(out, [[1,2,3], [4,5]])
end

function enz_enz(cat1::Symbol, cat2::Symbol)
  cat1sym = cat1
  cat2sym = cat2
  catcat = Symbol(cat1sym, sep_sym, cat2sym)
  out = oapply(enzXY, Dict(:bindXY=>bindunbind(cat1, cat2), :degXY=>degrade(cat1, cat2), :bindXYinact=>bindunbind(cat1, Symbol(cat2, :_inact)), :degXYinact=>degrade(cat1, Symbol(cat2, :_inact))), Dict(
    :X=>ob(cat1),
    :Xinact=>ob(Symbol(cat1sym,:_inact)),
    :Xdeg=>ob(Symbol(cat1sym,:_deg)),
    :Y=>ob(cat2),
    :Yinact=>ob(Symbol(cat2sym,:_inact)),
    :Ydeg=>ob(Symbol(cat2sym,:_deg)),
    :XY=>ob(catcat),
    :XYinact=>ob(Symbol(catcat,:_inact))))
  bundle_legs(out, [[1,2,3], [4,5,6]])
end

function enzyme_uwd(enzymes::Array{Symbol}, substrates::Array{Symbol})
  rel = RelationDiagram{Symbol}(0)

  chemicals = vcat(substrates, enzymes)

  subs = add_junctions!(rel, length(substrates), variable=substrates)
  enzs = add_junctions!(rel, length(enzymes), variable=enzymes)
  nsubs = length(subs)
  nenzs = length(enzs)

  catx = add_parts!(rel, :Box, nenzs, name=[Symbol("cat$i") for i in enzymes])
  add_parts!(rel, :Port, nenzs, junction=enzs, box=catx)

  for x in 1:nenzs
    for y in 1:nenzs
      if y != x
        catxy = add_part!(rel, :Box, name=Symbol("cat$(enzymes[x])cat$(enzymes[y])"))
        add_parts!(rel, :Port, 2, junction=[enzs[x], enzs[y]], box=catxy)
      end
    end
  end

  for x in 1:nenzs
    for y in 1:nsubs
      catxy = add_part!(rel, :Box, name=Symbol("cat$(enzymes[x])sub$(substrates[y])"))
      add_parts!(rel, :Port, 2, junction=[enzs[x], subs[y]], box=catxy)
    end
  end
  add_parts!(rel, :OuterPort, length(chemicals), outer_junction = vcat(subs, enzs))
  rel
end

function enzyme_generators(enzymes::Array{Symbol}, substrates::Array{Symbol})
  gens = Dict{Symbol, Any}()
  for e1 in enzymes
    for e2 in enzymes
      if e1 == e2
        gens[Symbol(:cat, e1)] = enz(e1)
      else
        gens[Symbol(:cat, e1, :cat, e2)] = enz_enz(e1, e2)
      end
    end
    for s in substrates
      gens[Symbol(:cat, e1, :sub, s)] = enz_sub(e1, s)
    end
  end
  gens
end
end

# ╔═╡ db1cf041-8f29-4433-a8b7-da5ea4828d52
Pkg.resolve()

# ╔═╡ aa5b099e-771f-42ff-a1b1-d22f0e615395


# ╔═╡ 178e764e-e239-4689-bb2f-4993b7755724
import .EnzymeReactions: ob, ode,
		   inactivate, bindunbind, degrade,
		   enzX, enzXY, enzXsubY,
		   enz, enz_enz, enz_sub,
		   enzyme_uwd, enzyme_generators

# ╔═╡ 563cf0a2-80e0-4bc2-8f6f-6a47cb2112af
@present ThAffinityNet(FreeSchema) begin
  (Rxn, Species)::Ob
  (Rate, Label)::Data
  binder::Hom(Rxn, Species)
	bindee::Hom(Rxn, Species)

  affinity::Attr(Rxn, Rate)

  family::Attr(Species, Label)
  form::Attr(Species, Label)
end;

# ╔═╡ 2d89b8e5-31a0-402c-b95c-87494a5a1317
begin
const AbstractAffinityNet = AbstractACSetType(ThAffinityNet)
const AffinityNet = ACSetType(ThAffinityNet){Real, Symbol}

function migrate_species!(an::AffinityNet, ln::LabelledReactionNet, index::Int64)
  labels = Symbol.(split(string(ln[index, :sname]), "_"))
  if length(labels) == 1
    push!(labels, :norm)
  end
  @assert (length(labels) == 2)  begin
    error("Label \"$(ln[index,:sname])\" is not in a valid format (family_form)")
  end
  add_part!(an, :Species, family=labels[1], form=labels[2])
end

function AffinityNet(ln::LabelledReactionNet{R, C}) where {R,C}
  # Collect associative/dissociative pairs
  assoc  = Dict{Int64, Array{Tuple{Int64, Int64, Int64},1}}()
  dissoc = Dict{Int64, Array{Tuple{Int64, Int64, Int64},1}}()
  for o in 1:nparts(ln, :T)
    outputs = ln[incident(ln, o, :ot), :os]
    inputs = ln[incident(ln, o, :it), :is]

    if length(outputs) == 2 && length(inputs) == 1
      if !(inputs[1] in keys(dissoc))
        dissoc[inputs[1]] = Array{Tuple{Int64, Int64, Int64},1}()
      end
      push!(dissoc[inputs[1]], outputs[1] < outputs[2] ? (outputs[1], outputs[2], o) :
                                                         (outputs[2], outputs[1], o))
    elseif length(inputs) == 2 && length(outputs) == 1
      if !(outputs[1] in keys(assoc))
        assoc[outputs[1]] = Array{Tuple{Int64, Int64, Int64},1}()
      end
      push!(assoc[outputs[1]], inputs[1] < inputs[2] ? (inputs[1], inputs[2], o) :
            (inputs[2], inputs[1], o))
    end
  end

  # Generate affinity net from pairs
  an = AffinityNet()
  s2i = Dict{Int64, Int64}()
  for prod in keys(assoc)
    if prod in keys(dissoc)
      for j in 1:length(assoc[prod])
        for i in 1:length(dissoc[prod])
          if assoc[prod][j][[1,2]] == dissoc[prod][i][[1,2]]
            ls, rs = ln[incident(ln, assoc[prod][j][3], :it), :is]
            ar, dr = (ln[assoc[prod][j][3], :rate], ln[dissoc[prod][i][3], :rate])
            if !(ls in keys(s2i))
              s2i[ls] = migrate_species!(an, ln, ls)
            end
            if !(rs in keys(s2i))
              s2i[rs] = migrate_species!(an, ln, rs)
            end
            add_part!(an, :Rxn, binder=s2i[ls], bindee=s2i[rs], affinity=ar/dr)
          end
        end
      end
    end
  end
  an
end

function propertygraph(an::AffinityNet;
                       prog::AbstractString="neato", graph_attrs::AbstractDict=Dict(),
                       node_attrs::AbstractDict=Dict(), edge_attrs::AbstractDict=Dict(:len=>"2.0"),
                       node_labels::Bool=true, edge_labels::Bool=false)

  shapes = ["box", "circle", "triangle", "diamond", "pentagon", "hexagon", "septagon", "octagon"]

  families = unique(an[:family])
  forms =	unique(an[:form])

  colors = string.(hex.(distinguishable_colors(length(forms), [RGB(1,1,1),RGB(0,0,0)], dropseed=true)))

  fam2shp = Dict(families[i]=>shapes[i] for i in 1:length(families))
  frm2clr = Dict(forms[i]=>colors[i] for i in 1:length(forms))

  rates = an[:affinity]
  minrate, maxrate = (minimum(rates), maximum(rates))
  rate_to_width(rate) = begin
    0.1 + (log(rate) - log(minrate)) * (2/(log(maxrate)-log(minrate)))
  end

	node_labeler(v) = begin
    Dict(:label=>"$(an[v,:family])$(an[v,:form]==:norm ? "" : Symbol("_",an[v,:form]))",
         :style=>"filled",
         :width=>"0.75", :height=>"0.75", :fixedsize=>"false",
         :shape=>fam2shp[an[v,:family]],
         :fillcolor=>"#$(frm2clr[an[v,:form]])")
  end

  edge_labeler(e) = begin
    return Dict(:color=>"black", :penwidth=>"$(round(rate_to_width(an[e,:affinity]), digits=1))")
  end

  g = Catlab.Graphs.Graph()
  migrate!(g, an, Dict(:V=>:Species, :E=>:Rxn), Dict(:tgt=>:bindee, :src=>:binder))

  Catlab.Graphs.PropertyGraph{Any}(g, node_labeler, edge_labeler;
                                  prog = prog,
                                  graph = merge!(Dict(:rankdir => "TB"), graph_attrs),
                                  node = merge!(Graphics.GraphvizGraphs.default_node_attrs(node_labels), node_attrs),
                                  edge = merge!(Dict(:arrowsize => "0.5"), edge_attrs),
                                  )
end

function Graphics.to_graphviz(an::AffinityNet; kwargs...)
  propertygraph(an; kwargs...) |> to_graphviz
end
end

# ╔═╡ 3779b846-e5ec-4239-a1d4-af2f8c2f10eb
begin
  using PlutoUI
  using LabelledArrays

  using DifferentialEquations
  using Plots
  plotly()

  display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>".75"));
  nothing
end

# ╔═╡ 93df89f0-8429-4fcc-bd01-6982417f5134
begin
  # Initial Concentrations
  K = :K=>33000;
  S = :S=>33000;
  L = :L=>33000;
  Kinact = :K_inact=>0;
  Sinact = :S_inact=>0;
  Linact = :L_inact=>0;
  E = :E=>700000;
  G = :G=>714000;
  Spike = :Spike=>66000;

  cats = [:K, :S, :L]
  subs = [:E, :G, :Spike]

  # Parameter Rates (units of pM and min)
  rxns = Dict(
    "K" => [inactivate(K, 7.494e-10)
           bindunbind(K, K, 7.814e-4, 3.867e-3)
           degrade(K, K, 2.265e-1)
           bindunbind(K, Kinact, 7.814e-4, 3.867e-3)
           degrade(K, Kinact, 2.265e-1)],
    "S" => [inactivate(S, 1.906e-2)
           bindunbind(S, S, 3.534e-7, 1.688e2)
           degrade(S, S, 1.433e1)
           bindunbind(S, Sinact, 3.534e-7, 1.688e2)
           degrade(S, Sinact, 1.433e1)],
    "L" => [inactivate(L, 7.810e-3)
           bindunbind(L, L, 1.000e-11, 7.440e3)
           degrade(L, L, 2.670e2)
           bindunbind(L, Linact, 1.000e-11, 7.440e3)
           degrade(L, Linact, 2.670e2)],
    "K↦E" => [bindunbind(K, E, 9.668e-6,1.000e-2)
            degrade(K, E, 1.728e0)],
    "K↦G" => [bindunbind(K, G, 2.764e-6, 8.780e-1)
            degrade(K, G, 1.502e0)],
    "S↦E" => [bindunbind(S, E, 4.197e-7, 1.06e-3)
            degrade(S, E, 1.384e4)],
    "S↦G" => [bindunbind(S, G, 5.152e-8, 3.894e-3)
            degrade(S, G, 8.755e-1)],
    "L↦E" => [bindunbind(L, E, 1.977e-8, 1.000e-2)
            degrade(L, E, 1.066e2)],
    "L↦G" => [bindunbind(L, G, 3.394e-8, 2.365e1)
            degrade(L, G, 4.352e0)],
    "K↦S" => [bindunbind(K, S, 8.822e-4, 4.114e5)
            degrade(K, S, 9.000e-10)
            bindunbind(K, Sinact, 8.822e-4, 4.114e5)
            degrade(K, Sinact, 9.000e-10)],
    "K↦L" => [bindunbind(K, L, 1.756e-4, 3.729e4)
            degrade(K, L, 6.505e6)
            bindunbind(K, Linact, 1.756e-4, 3.729e4)
            degrade(K, Linact, 6.505e6)],
    "S↦K" => [bindunbind(S, K, 3.679e-4, 1.562e3)
            degrade(S, K, 4.410e2)
            bindunbind(S, Kinact, 3.679e-4, 1.562e3)
            degrade(S, Kinact, 4.410e2)],
    "S↦L" => [bindunbind(S, L, 1.000e-3, 5.000e2)
            degrade(S, L, 1.000e-7)
            bindunbind(S, Linact, 1.000e-3, 5.000e2)
            degrade(S, Linact, 1.000e-7)],
    "L↦K" => [bindunbind(L, K, 1.000e-3, 4.118e3)
            degrade(L, K, 3.234e1)
            bindunbind(L, Kinact, 1.000e-3, 4.118e3)
            degrade(L, Kinact, 3.234e1)],
    "L↦S" => [bindunbind(L, S, 1.056e-12, 5.000e2)
            degrade(L, S, 5.000e-1)
            bindunbind(L, Sinact, 1.056e-12, 5.000e2)
            degrade(L, Sinact, 5.000e-1)]
  );
  rxns = Dict(Symbol(k)=>v for (k,v) in rxns)
  # define labels to reaction network mappings
  functor(x) = oapply(x, Dict(
    :catK=>enz(rxns, K),
    :catS=>enz(rxns, S),
    :catL=>enz(rxns, L),
    :catKcatS=>enz_enz(rxns, K,S),
    :catKcatL=>enz_enz(rxns, K,L),
    :catScatK=>enz_enz(rxns, S,K),
    :catScatL=>enz_enz(rxns, S,L),
    :catLcatK=>enz_enz(rxns, L,K),
    :catLcatS=>enz_enz(rxns, L,S),
    :catKsubE=>enz_sub(rxns, K,E),
    :catSsubE=>enz_sub(rxns, S,E),
    :catLsubE=>enz_sub(rxns, L,E),
    :catKsubG=>enz_sub(rxns, K,G),
    :catSsubG=>enz_sub(rxns, S,G),
    :catLsubG=>enz_sub(rxns, L,G)));
  lfunctor(x) = oapply(x, enzyme_generators([:K,:S,:L],[:G,:E,:Spike]))
  def_model = apex(functor(enzyme_uwd([:K,:S,:L],[:G,:E])))
  def_rates = rates(def_model)
  def_concs = Dict(c=>concentrations(def_model)[c] for c in snames(def_model))
  def_concs[:Spike] = Spike[2]
  def_concs[:Spike_deg] = 0
  def_concs[Symbol("K↦Spike")] = 0
  def_concs[Symbol("S↦Spike")] = 0
  def_concs[Symbol("L↦Spike")] = 0
  nothing
end

# ╔═╡ fe9b889d-79c2-493b-9426-e33e6820cd90
md""" Upload a rate data file to use:  $(@bind user_csv FilePicker()) """

# ╔═╡ 950d3b4e-f957-45b6-aa80-e3dfc765aad0
md"""Use default rates (override uploaded rates)? $(@bind defbox CheckBox(false))"""

# ╔═╡ 50334069-a50c-467c-94ae-63b9b2264a18
begin
	local k_rates, status
	try
		defbox && error("default rate override")
		k_rates = UInt8.(user_csv["data"]) |> IOBuffer |> JSON.parse
		status = md"""Using rates from $(user_csv["name"])"""
	catch e
		k_rates = def_rates
		status = md"""No file selected, reverting to default rates"""
	end
	m_rates = Dict(Symbol(k)=>k_rates[k] for k in keys(k_rates))
	valid_enz_sub = Set{Symbol}()
	for e in cats
		for s in vcat(cats, subs)
      if Symbol("bind_$e↦$s") ∈ keys(m_rates)
				push!(valid_enz_sub, e)
				push!(valid_enz_sub, s)
			end
		end
	end
	status
end

# ╔═╡ ddc141ba-d2e8-4ac4-8bc3-12fb1bb9fd4d
md"""
#### Enzyme Diagram
Enzymes

"""

# ╔═╡ 4ad16c5c-73bc-4e42-9bfc-aea73a6bfbfe
begin
#TODO
#Find a more compact way to represent these checkboxes while still dynamically #determining which to display

	local boxes = Dict(:K=>md"K $(@bind kbox CheckBox(:K ∈ valid_enz_sub))",
				 :S=>md"S $(@bind sbox CheckBox(:S ∈ valid_enz_sub))",
				 :L=>md"L $(@bind Lbox CheckBox(:L ∈ valid_enz_sub && false))")
md"""
$(:K ∈ valid_enz_sub ? boxes[:K] : md"")
$(:S ∈ valid_enz_sub ? boxes[:S] : md"")
$(:L ∈ valid_enz_sub ? boxes[:L] : md"")
"""
end

# ╔═╡ e89794b1-5bcd-4b6c-9cb2-77deca569c2e
md"""Substrates"""

# ╔═╡ dcdb88ef-f04f-4ee8-87cc-bb26f396f064
begin

	local boxes = Dict(:G=>md"G $(@bind gbox CheckBox(:G ∈ valid_enz_sub))",
				 :E=>md"E $(@bind ebox CheckBox(:E ∈ valid_enz_sub && false))",
				 :Spike=>md"Spike $(@bind pbox CheckBox(:Spike ∈ valid_enz_sub))")
md"""
$(:G ∈ valid_enz_sub ? boxes[:G] : md"")
$(:E ∈ valid_enz_sub ? boxes[:E] : md"")
$(:Spike ∈ valid_enz_sub ? boxes[:Spike] : md"")
"""
end

# ╔═╡ d9f5de8a-f3a2-41c9-9f3c-a0c8347368a4
begin
inp = [kbox, sbox, Lbox];
inpSymb = [:K, :S, :L];
inpF = Symbol[];
	for i = 1:3
		if inp[i]
			push!(inpF,inpSymb[i]);
		end

	end

outp = [gbox, ebox, pbox];
outSymb = [:G, :E, :Spike];
outF = Symbol[];
	for i = 1:3
		if outp[i]
			push!(outF,outSymb[i]);
		end

	end

end


# ╔═╡ e6589d31-dce7-42c3-b494-db03fe561ae9
  uwd = enzyme_uwd(inpF, outF);

# ╔═╡ 7dbe9349-8b9e-4ac2-b4bf-b59f58a10ebc
begin
  display_uwd(uwd)
end

# ╔═╡ cf9e03db-42b7-41f6-80ce-4b12ddb93211
begin 
model = uwd |> lfunctor |> apex;
  r = join(["'$k': $(m_rates[k])" for k in tnames(model)], ", ");

  r2 = join(["'$k': $(def_concs[k])" for k in snames(model)], ", ");
  #== Basic Slider Template
<li><input
  class='slider'
  type='range'
  min='-10'
  step='0.01'
  max='10'
  value='0'
  oninput='this.nextElementSibling.value=(10**this.value).toExponential(2)'
  ></input>
<output>1</output></li>
==#
c = vcat(["$(m_rates[k])" for k in tnames(model)], ["$(def_concs[k])" for k in snames(model)]);

form_vals = HTML("""
<bond def = "c">
<form>
  <div id ="myDIV" style='height:500px;overflow:scroll'>
  <h4>Adjust Rates</h4>
  <table>
  </table>
  </div>
 <div id ="concDiv" style='height:500px;overflow:scroll'>
		<h4>Adjust Concentrations</h4>
		<table id = "table2">
		</table>
	</div>
 <div id ="bindingAff" style='height:500px;overflow:scroll'>
		<h4>Binding Affinities</h4>
		<table id = "bindingAffTable">
		</table>
	</div>
</form>
</bond>
<style>
#myDIV {
  width: 100%;
  padding: 50px 0;
  text-align: center;
 background-color: #F9F6E5;
  margin-top: 20px;
}
#concDiv {
  width: 100%;
  padding: 50px 0;
  text-align: center;
  background-color:  #F9F6E5;
  margin-top: 20px;
  display: none
}
#bindingAff {
  width: 100%;
  padding: 50px 0;
  text-align: center;
  background-color:  #F9F6E5;
  margin-top: 20px;
  display: none
}
.button {
  background-color: #003057;
  border: none;
  color: white;
  padding: 15px 32px;
  text-align: center;
  text-decoration: none;
  display: inline-block;
  font-size: 16px;
  margin: 4px 2px;
  cursor: pointer;
}
.button:focus {
    background-color:#BFB37C;
}
.slider {
  -webkit-appearance: none;
  width: 100%;
  height: 7px;
  border-radius: 5px;
  background: #d3d3d3;
  outline: none;
  opacity: 0.7;
  -webkit-transition: .2s;
  transition: opacity .2s;
}
.slider::-webkit-slider-thumb {
  -webkit-appearance: none;
  appearance: none;
  width: 20px;
  height: 20px;
  border-radius: 50%;
  background:  #BFB37C;
  cursor: pointer;
}
.slider::-moz-range-thumb {
  width: 20px;
  height: 20px;
  border-radius: 50%;
  background: #003057;
  cursor: pointer;
}
</style>
<script>
//`currentScript` is the current script tag - we use it to select elements//
function hideShowRate() {
	var x = document.getElementById("myDIV");
	var y = document.getElementById("concDiv");
	var z = document.getElementById("bindingAff");
	z.style.display = "none";
	x.style.display = "block";
	y.style.display = "none"
}
function hideShowConc() {
	var x = document.getElementById("concDiv");
	var y = document.getElementById("myDIV");
	var z = document.getElementById("bindingAff");
	x.style.display = "block";
	y.style.display = "none"
	z.style.display = "none"
}
function hideShowBind() {
	var x = document.getElementById("bindingAff");
	var y = document.getElementById("myDIV");
	var z = document.getElementById("concDiv");
	x.style.display = "block";
	y.style.display = "none"
	z.style.display = "none"
	generateBindingTable()
}
function returnStrings(s,arr){
		let finalArr = []
		for (let item in arr){
			if(arr[item].includes(s) == true){
				finalArr.push(arr[item])
				}
			}
		return finalArr
		}
function generateBindingTable(){
		let keyNames = Object.keys(rates)
		let objVals = Object.values(rates)
		let bindArr = []
		let unbindArr = []
		for (let i in keyNames){
		   let nm = keyNames[i]
			if(nm.substring(0,4) == 'bind'){
				bindArr.push(nm)
		}else if(nm.substring(0,6) == 'unbind'){
				unbindArr.push(nm)
			}
			}
		let kdArr = []; //array of binding affinities
		let kdNames = []; // array of binding affinity names
			//grab values from sliders
	    let x = form.getElementsByClassName('sliderval');
		let xarr = Array.from(x, (v,_)=>{return v.value})
		//Loop over bind and unbind arrays to do calculations
		for(let i in keyNames){
			if(keyNames[i].substring(0,4) == 'bind'){
			let bind = keyNames[i]
			for(let j in keyNames){
		    if(keyNames[j].substring(0,6) == 'unbind'){
				let unbind = keyNames[j]
					if(bind.substring(5,bind.length) == unbind.substring(7,unbind.length)){
								kdArr.push(xarr[j]/xarr[i])
						kdNames.push('Kd_'+bind.substring(5,bind.length))
		//Create name of deg
		let nm = bind.substring(4,bind.length)
		let deg_name = 'deg' + nm
			for(let k in keyNames){
				if(keyNames[k]== deg_name){
						let kcat = keyNames[k]
						let unbindval = parseFloat(xarr[j]);
						let bindval = parseFloat(xarr[i])
						let catval = parseFloat(xarr[k])
						let k_m = (catval+unbindval)/bindval;
						kdArr.push(k_m)
						kdNames.push('Km_'+bind.substring(5,bind.length))
						kdArr.push(catval/k_m)
						kdNames.push('Cat_efficieny_'+bind.substring(5,bind.length))
							}
						}}}
					}
				}
			}
		const list3 = form.querySelector('#bindingAffTable')
		while(list3.rows.length > 0) {
				  list3.deleteRow(0);
				}
		var header = document.createElement('thead')
 		var toprow = document.createElement('tr');
		var toprowItem = document.createElement('th');
		var toprowItem2 = document.createElement('th');
		var toprowItem3 = document.createElement('th');
		var toprowItem4 = document.createElement('th');
		toprowItem.innerText = 'Name';
		toprowItem2.innerText = 'Kd';
		toprowItem3.innerText = 'Km';
		toprowItem4.innerText = 'Catalytic Efficiency';
		header.appendChild(toprow)
		//toprow.appendChild(header);
		toprow.appendChild(toprowItem);
		toprow.appendChild(toprowItem2);
		toprow.appendChild(toprowItem3);
		toprow.appendChild(toprowItem4);
		list3.appendChild(toprow);
		for (let i = 0; i < kdNames.length/3; i++){
		  var item = document.createElement('tr');
		  var label = document.createElement('td');
		  label.innerText = kdNames[i*3].substring(3)
		  var label2 = document.createElement('td');
		  var label3 = document.createElement('td');
		  var label4 = document.createElement('td');
		  label2.innerText =  kdArr[i*3].toExponential(2)
		label3.innerText = kdArr[i*3+1].toExponential(2) // 'placeholder'
		label4.innerText = kdArr[i*3+2].toExponential(2) // 'placeholder'
		  item.appendChild(header)
		  item.appendChild(label)
          item.appendChild(label2)
		  item.appendChild(label3)
		  item.appendChild(label4)
		  list3.appendChild(item)
		}
		}
const form = currentScript.parentElement.querySelector('form')
const list = form.querySelector('table')
var rates = {$r};
var conc = {$r2};
for ( var r in rates ){
  console.log(r);
  var item = document.createElement('tr');
  var label = document.createElement('th');
  label.innerText = r
  var slider_box = document.createElement('td');
  var slider_val_box = document.createElement('td');
  var slider = document.createElement('input');
  var slider_val = document.createElement('input');
  slider_val.setAttribute('type', 'text');
  slider_val.setAttribute('class', 'sliderval');
  slider_val.value = rates[r].toExponential(2);
  slider.setAttribute('type', 'range');
  slider.setAttribute('class', 'slider');
  slider.setAttribute('min', '-9.99');
  slider.setAttribute('max', '10');
  slider.setAttribute('value', '0.0');
  slider.setAttribute('step', '0.01');
  slider.setAttribute('oninput', `this.parentElement.nextElementSibling.children[0].value=((10**this.value)*\${rates[r]}).toExponential(2)`);
  slider_box.appendChild(slider)
  slider_val_box.appendChild(slider_val)
  item.appendChild(label)
  item.appendChild(slider_box)
  item.appendChild(slider_val_box)
  list.appendChild(item)
}
const list2 = form.querySelector('#table2')
for ( var r in conc ){
  console.log(r);
  var item = document.createElement('tr');
  var label = document.createElement('th');
  label.innerText = r
  var slider_box = document.createElement('td');
  var slider_val_box = document.createElement('td');
  var slider = document.createElement('input');
  var slider_val = document.createElement('input');
  slider_val.setAttribute('type', 'text');
  slider_val.setAttribute('class', 'sliderval');
  slider_val.value = conc[r].toExponential(2);
  slider.setAttribute('type', 'range');
  slider.setAttribute('class', 'slider');
  slider.setAttribute('min', '-9.99');
  slider.setAttribute('max', '10');
  slider.setAttribute('value', '0.0');
  slider.setAttribute('step', '0.01');
  console.log('concentration')
  console.log(conc[r]);
  let test = conc[r];
  if(test == 0){
		test = .000001
		}
  slider.setAttribute('oninput', `this.parentElement.nextElementSibling.children[0].value=((10**this.value)*\${test}).toExponential(2)`);
  slider_box.appendChild(slider)
  slider_val_box.appendChild(slider_val)
  item.appendChild(label)
  item.appendChild(slider_box)
  item.appendChild(slider_val_box)
  list2.appendChild(item)
}
var x = form.getElementsByClassName('sliderval');
function onsubmit(){
  // We send the value back to Julia //
  form.value = Array.from(x, (v,_)=>{return v.value})
  form.dispatchEvent(new CustomEvent('input'))
  console.log(form.value)
}
var b = document.createElement('input');
b.setAttribute('type', 'button');
b.setAttribute('class', 'button');
b.value = 'Update Plot';
b.addEventListener('click', function() {onsubmit();console.log('hello from button')})
form.appendChild(b)
onsubmit()
let sp2 = document.getElementById("myDIV")
var d = document.createElement('input')
d.setAttribute('type','button')
d.setAttribute('class','button')
d.value = 'Rates'
d.addEventListener('click', function() {hideShowRate();})
form.insertBefore(d,sp2)
var e = document.createElement('input')
e.setAttribute('type','button')
e.setAttribute('class','button')
e.value = 'Concentrations'
e.addEventListener('click', function() {hideShowConc();})
form.insertBefore(e,sp2)
var f = document.createElement('input')
f.setAttribute('type','button')
f.setAttribute('class','button')
f.value = 'Binding Affinities'
f.addEventListener('click', function() {hideShowBind();})
form.insertBefore(f,sp2)
</script>
""");
end

# ╔═╡ 1f27f695-9d83-4744-b6f0-e3280b746691
length(c)

# ╔═╡ bc90290e-81ba-4963-8a1f-5be15c0f3940
length(tnames(model))

# ╔═╡ 9e847b90-54ae-4afb-bf72-1b9adeb97abc
length(snames(model))

# ╔═╡ 066b7505-e21b-467e-86c1-cea1ff80246e
begin
cur_rate = Dict(tnames(model)[i]=>parse(Float64, c[i]) for i in 1:length(tnames(model)))
  # cur_conc = concentrations(model)
  cur_conc_prep =  Dict(snames(model)[i]=>parse(Float64, c[i+length(tnames(model))]) for i in 1:length(snames(model)))

	keysa = keys(cur_conc_prep)
	keyArr = Symbol[]
	keyArrStr= String[]
	valArr = Float64[]
	for key in keys(cur_conc_prep)
		push!(keyArr,key)
		push!(valArr,cur_conc_prep[key])
		push!(keyArrStr,String(key))
	end

  cur_conc = @LArray valArr Tuple(keyArr)
  vf = vectorfield(model);
sol = solve(ODEProblem(vf, cur_conc, (0.0,120.0),cur_rate));
  nothing
end

# ╔═╡ 3f8db202-ac50-462d-b96d-ba629ca43325
md"""### Iterate over multiple species concentrations $(@bind multCheck CheckBox(false))
#### Select species to iterate: 
$(@bind veg Select(keyArrStr))

Starting Concentration: $(@bind initConc TextField()) \
End Concentration: $(@bind endConc TextField()) \
number of steps: $(@bind stepsConc NumberField(1:10, default=5)) \

Choose spacing type: $(@bind spacing Select(["Logarithmic","Linear"])) \
Display combined concentrations excluding degraded? $(@bind combineCheckDeg CheckBox(false)) \
Use constant scale for axis? $(@bind scaleCheck CheckBox(false))

"""

# ╔═╡ cb45ae16-8893-4b80-bae0-fcef646a56b2
begin
if multCheck
keyArr3 = ["A","B"];
graphKeys3 = []
graphKeyVals3 = HTML("""
<bond def = "graphKeys3">
<form>

  <div id ="myDIV3" style='height:300px;overflow:scroll'>

  <h4>Select Species to plot</h4>
			<p>This plot will display how selected species change over iterations. Limit of 6 species. Iteration is indicated by color and species is indicated by marker shape</p>
			
  <select  id="graphConc" multiple = "multiple" size = "10">

  </select>



  </div>



</form>
</bond>
<style>

#myDIV2 {
  width: 100%;
  padding: 10px 0;
  text-align: center;
  background-color: #F9F6E5;
  margin-top: 10px;
}


</style>



<script>
//`currentScript` is the current script tag - we use it to select elements//


const form = currentScript.parentElement.querySelector('form')
const list =  form.querySelector('select')
console.log(list)
var concList = $keyArrStr;
console.log(concList)

for ( var key in concList ){
  console.log(concList[key]);
  //var item = document.createElement('tr');
  //var label = document.createElement('th');
  var item = document.createElement('option');
  item.value = concList[key];
  item.label = concList[key];
  list.appendChild(item)

}





var x = form.getElementsByClassName('selector');
function onsubmit(){
		console.log('onsubmit called')
  // We send the value back to Julia //
		console.log(x)


    var selected = [];
    for (var option of document.getElementById('graphConc').options)
    {
        if (option.selected) {
            selected.push(option.value);
        }
    }
 form.value = selected
  form.dispatchEvent(new CustomEvent('input'))

  console.log(form.value)
    console.log(selected)
}
var b = document.createElement('input');
b.setAttribute('type', 'button');
b.setAttribute('class', 'button');
b.value = 'Plot species';
b.addEventListener('click', function() {onsubmit();console.log('hello from button')})
form.appendChild(b)
onsubmit()




</script>
""");
	else
		nothing;
	end
end

# ╔═╡ 12866252-a5c6-43d0-92f1-d52df5a2d949
md"""Display total concentrations? $(@bind combineCheck CheckBox(false)) 

"""

# ╔═╡ d80f94c4-03d2-4aac-90f5-9415405b4412
begin

keyArr2 = ["A","B"];
graphKeys = []
graphKeyVals = HTML("""
<bond def = "graphKeys">
<form>

  <div id ="myDIV2" style='height:300px;overflow:scroll'>

  <h4>Select Variables to Export</h4>
  <select  id="graphConc" multiple = "multiple" size = "10">

  </select>



  </div>



</form>
</bond>
<style>

#myDIV2 {
  width: 100%;
  padding: 10px 0;
  text-align: center;
  background-color: #F9F6E5;
  margin-top: 10px;
}


</style>



<script>
//`currentScript` is the current script tag - we use it to select elements//


const form = currentScript.parentElement.querySelector('form')
const list =  form.querySelector('select')
console.log(list)
var concList = $keyArrStr;
console.log(concList)

for ( var key in concList ){
  console.log(concList[key]);
  //var item = document.createElement('tr');
  //var label = document.createElement('th');
  var item = document.createElement('option');
  item.value = concList[key];
  item.label = concList[key];
  list.appendChild(item)

}





var x = form.getElementsByClassName('selector');
function onsubmit(){
		console.log('onsubmit called')
  // We send the value back to Julia //
		console.log(x)


    var selected = [];
    for (var option of document.getElementById('graphConc').options)
    {
        if (option.selected) {
            selected.push(option.value);
        }
    }
 form.value = selected
  form.dispatchEvent(new CustomEvent('input'))

  console.log(form.value)
    console.log(selected)
}
var b = document.createElement('input');
b.setAttribute('type', 'button');
b.setAttribute('class', 'button');
b.value = 'Update variable list';
b.addEventListener('click', function() {onsubmit();console.log('hello from button')})
form.appendChild(b)
onsubmit()




</script>
""");

end

# ╔═╡ 043f7a23-3b59-4e34-a8d3-9853cc66c228
md"""Filter Inact? $(@bind filterInact CheckBox(false)) \
Filter degraded? $(@bind filterDeg CheckBox(false))
"""


# ╔═╡ 1596bc9f-f7e4-4d3d-9978-9da4eecbaede
begin
	filter_list = [];
	
if filterInact == true || filterDeg == true
	if filterInact == true
			append!(filter_list,["inact"]);
	end
	if filterDeg == true
			append!(filter_list,["deg"]);
	end
		
end
	nothing;
end

# ╔═╡ 675d0bb0-4601-4f4e-bc7d-5d5fb2d70b18
md"""Export only selected variables? $(@bind importCheck CheckBox(false))"""

# ╔═╡ ad8edd69-c164-4221-bdee-e7c9381ffcab
begin

graphKeySymb = Symbol[]
for item in graphKeys
		push!(graphKeySymb,Symbol(item))

end
end

# ╔═╡ 56afefe8-4452-4b2a-8a3b-e493ee1dd6c6
begin
	function formatSymbArr(arr)
		#Change string array to symbol
		B = Array{Symbol}(undef, length(arr))
		for i in 1:length(arr)
			 B[i] = Symbol(arr[i])
		end
		return B
	end


	function formatStrArr(arr)
		#Change symbol array to string
		B = Array{String}(undef, length(arr))
		for i in 1:length(arr)
			 B[i] = string(arr[i])
		end
		return B

	end

	function combineConc(df)
		#Build a list of the different chemicals. Basically take the first letter of each and the make the list unique so that we can combine everything easily.
		cols = names(df)
		# println(cols)
		chems = Array{String}(undef, 0)
		for i in 2:length(cols)

			# letter = cols[i][1]
			letter = split(cols[i],"↦")[1]
			println(split(cols[i],"↦")[1])
			append!(chems,[letter])
		end
		unique!(chems)
		println(chems)
		dff = DataFrame()
		dff.timestamp = df[!,"timestamp"]
		for i in 1:length(chems)
			arr = zeros(length(dff.timestamp))
			# println("new chem")
			for j in 2:length(cols)
				letter2 = split(cols[j],"↦")[1]
				# if chems[i] == cols[j][1]
				println(letter2)
				if chems[i] == letter2
					# println("arr before: ",arr[10] )
					# println("Value added: ", df[10,cols[j]])
					arr = arr .+ df[!,cols[j]]
					# println("arr after: ",arr[10] )

				end
			end

			colname = string(chems[i])

			dff[!,colname] = arr
		end

		return dff
	end
	
	
		function combineConcNoDeg(df)
		#Build a list of the different chemicals. Basically take the first letter of each and the make the list unique so that we can combine everything easily.
		cols = names(df)
		# println(cols)
		chems = Array{String}(undef, 0)
		for i in 2:length(cols)

			# letter = cols[i][1]
			letter = split(cols[i],"↦")[1]
			# println(split(cols[i],"↦")[1])
			if occursin("deg",letter) == false
				append!(chems,[letter])
			end
		end
		unique!(chems)
		# println(chems)
		dff = DataFrame()
		dff.timestamp = df[!,"timestamp"]
		for i in 1:length(chems)
			arr = zeros(length(dff.timestamp))
			# println("new chem")
			for j in 2:length(cols)
				letter2 = split(cols[j],"↦")[1]

				if chems[i] == letter2
					arr = arr .+ df[!,cols[j]]

				end
			end

			colname = string(chems[i])

			dff[!,colname] = arr
		end

		return dff
	end
	
	
	function logrange(x1, x2, n) 
		B = Float64[] 
		startLog = log10(parse(Float64,x1))
		endLog = log10(parse(Float64,x2))
		
		logvals = range(startLog,endLog,length = n)
		for i in logvals
			append!(B,10.0^i)	
		end
		return B
		# return (10^y for y in range(log10(parse(Float64,x1)), log10(parse(Float64,x2)), length=n))
		
	end
	
	function createMultiArr(cur_conc,conList,veg)
		multiConcArr = []
		
		for i in 1:length(conList)
			tempArr = copy(cur_conc)
			println(tempArr)
			item = conList[i]
			println(item)
			tempArr[Symbol(veg)] = item
			println(tempArr)
			append!(multiConcArr,[tempArr])
		end
		return multiConcArr
	end
	
	function createPlot(sol, model,current_conc,n)
		tsteps = sol.t
		labels = isempty(graphKeySymb) ? snames(model) : graphKeySymb
		name = string("Plot ",string(n)," init conc: ", string(current_conc))
	  plot(tsteps, [[sol(t)[l]/1e3 for t in tsteps] for l in labels], labels=hcat(String.(labels)...),title = name, linewidth=3, xlabel="Minutes", ylabel="Solution Concentration (nM)")
		
	end
	
	function createPlotDF(sol,model,current_conc,n)
			graphKeySymb2 = Symbol[]
			labels2 = isempty(graphKeySymb2) ? snames(model) : graphKeySymb2
			
			S_labels2 = formatStrArr(labels2)
			prepend!(S_labels2,["timestamp"])
			dfs2 = DataFrame(sol)
			rename!(dfs2,S_labels2)
			# names(dfs)[1][1]
			# dff2 = combineConc(dfs2)
			dff2 = combineConcNoDeg(dfs2)
			labels_new2 = names(dff2)[2:end]
			name = string("Plot ",string(n),", init conc: ", string(current_conc))
			timesteps2 = dff2[!,"timestamp"]
			dff2[!,2:end]
			data2 = Matrix(dff2[!,2:end])/1e3
			return plot(timesteps2,data2, label = reshape(labels_new2, (1,length(labels_new2))), title = name, linewidth = 3, xlabel = "Minutes",ylabel = "Solution Concentration (nM)")
		
	end
	
	function generatePlotDF2(df_plot,timestamps,numSpecs,nplots)
		clrs = ["red" "blue" "green" "purple" "orange" "pink"]
		mkrs = [:circle :cross :rect :diamond :utriangle :dtriangle ]
		
		#Format markers and colors correctly
		clr_array = []
		mkrs_array = []
		for i in 1:nplots
			for j in 1:numSpecs
				append!(clr_array,[clrs[i]])
				append!(mkrs_array,[mkrs[j]])
			end
		end
		
		clr_array = reshape(clr_array, (1,length(clr_array)))
		mkrs_array = reshape(mkrs_array, (1,length(mkrs_array)))
		# mkrs = mkrs[1:numSpecs]
		# println(mkrs)
		mat_data = Matrix(df_plot)
		plotnames = names(df_plot)
		plotlabel = reshape(plotnames, (1,length(plotnames)))
		return plot(timestamps,mat_data, label = plotlabel, linecolor = clr_array, markershape = mkrs_array, markercolor = clr_array, linewidth = 3, xlabel = "Minutes",ylabel = "Solution Concentration (nM)")
	end
	
end



# ╔═╡ f7632708-7a0f-4c18-9f0b-44fb49aaeaa0
begin
if multCheck
			if spacing == "Linear"
		conList = range(parse(Float64,initConc), parse(Float64,endConc), length = stepsConc)
		
	elseif spacing == "Logarithmic"
		conList = logrange(initConc,endConc,stepsConc);
		# println(conList)
	else
		conList = ["Empty"]
	end
	multiConcArr = createMultiArr(cur_conc,conList,veg)
	sol_arr = []
	plotArr = []
	plotNumArr = []
	
	for item in multiConcArr
		tempsol = solve(ODEProblem(vf, item, (0.0,120.0),cur_rate), saveat = 1.0);
		append!(sol_arr,[tempsol])	
	end
		
	numPlots = length(multiConcArr)
	for p in 1:length(sol_arr)
		tsol = sol_arr[p]
		concName = conList[p]
		if combineCheckDeg == false
			append!(plotArr,[createPlot(tsol,model,concName,p)])
		else
			append!(plotArr,[createPlotDF(tsol,model,concName,p)])	
		end
		
	end
	
	# labels_mult = isempty(graphKeySymb) ? snames(model) : graphKeySymb
	labels_mult = collect(keys(multiConcArr[1]))
	for i in 1:length(conList)
		append!(plotNumArr,[string("Plot ",string(i))])
	end
		if scaleCheck == false

			plot(plotArr..., size = (600, 400*length(multiConcArr)),layout = (length(multiConcArr),1), legend = false, bottom_margin = 10mm)
		else
			maxlim = 1.0
			for j in 1:length(sol_arr)
				
				if maximum(sol_arr[j])>maxlim
					maxlim = maximum(sol_arr[j])
				end
			end
				plot(plotArr..., size = (600, 400*length(multiConcArr)),layout = (length(multiConcArr),1), legend = false, bottom_margin = 10mm, ylims = (0,maxlim/1e3))	
		end
	else
		plotNumArr = ["Empty"];
		nothing;
		
end
	
end

# ╔═╡ 99758072-3774-4bdc-9730-f0efff0e9330
begin 
	df_array = [];
if multCheck
		
		#Convert solutions to dataframe
		str_labels_mult = formatStrArr(labels_mult)
		numSpecs = length(graphKeys3)
		nplots = length(sol_arr)
		prepend!(str_labels_mult,["timestamp"])
		timestamps = []
	for i in 1:length(sol_arr)
			
		tempdf = DataFrame(sol_arr[i])
		rename!(tempdf,str_labels_mult)
		timestamps = tempdf[!, "timestamp"]
		append!(df_array,[select(tempdf,graphKeys3)])
	end
	df_plot = df_array[1]
	# timestamps = sol_arr[1][!,:timestamp]
	for i in 2:length(sol_arr)
		df_plot = hcat(df_plot,df_array[i], makeunique=true)
	end
		# println(df_plot)
		
		# mat_data = Matrix(df_plot)
		# plotnames = names(df_plot)
		# plotlabel = reshape(plotnames, (1,length(plotnames)))
		# plot(timestamps,mat_data, label = plotLabel, linewidth = 3, xlabel = "Minutes",ylabel = "Solution Concentration (nM)")
	else
		nothing;
	end
end

# ╔═╡ 95cfdfe4-7fd8-4155-8a29-5a9711ab0f10
begin 
	if multCheck == true && length(df_array)>0 
generatePlotDF2(df_plot,timestamps,numSpecs,nplots)
	else
		nothing;
	end
end

# ╔═╡ f0eb5d23-e4c5-4189-a4aa-adf3219227cf
md"""#### Choose export mode $(@bind exportMode Select(["Export all plots", "Export specific plot"]))
##### Select plot for export: 
$(@bind conListSel Select(formatStrArr(plotNumArr)))
"""

# ╔═╡ b9b223de-0ff8-436e-9a4c-e056e1c3a412
begin
	if multCheck
	if exportMode == "Export specific plot"
		ind = parse(Int8,conListSel[end])
		println("EXporting specific Plot")
		sol_sel = sol_arr[ind][1]
		df_sel = DataFrame(sol_arr[ind])
		nms_sel = collect(keys(sol_sel))
		prepend!(nms_sel,[:timestamp])
		rename!(df_sel,nms_sel)
		nmsStr_sel = formatStrArr(nms_sel)
		filtered_names_sel = []	
		filtered_names_inx_sel = []

		if length(filter_list) > 0
			for i in 1:length(filter_list)
				name = filter_list[i]
					
					for j in 1:length(nmsStr_sel)

						if occursin(name,nmsStr_sel[j]) == true
							append!(filtered_names_inx_Sel,[j])
						end
					end
			end
		end
		sort!(unique!(filtered_names_inx_sel))
		deleteat!(nmsStr_sel,filtered_names_inx_sel)

		finKeys_sel = formatSymbArr(nmsStr_sel)

		dfFin_sel = select(df_sel,finKeys_sel)
		# CSV.write("sim_res_sel.csv", dfFin_sel)
		XLSX.writetable("sim_results.xlsx", collect(DataFrames.eachcol(dfFin_sel)), DataFrames.names(dfFin_sel);overwrite = true)

		
	elseif exportMode == "Export all plots"
		println("Export all plots")
		# Use the xlsx package to create an excel sheet with multiple tabs. Then can export dataframes
		#First convert solution into dataframes
		df_list = []

		for sltn in 1:length(plotNumArr)
			println(sltn)
			temp_sol = sol_arr[sltn]
			# println(temp_sol)
			df_temp = DataFrame(temp_sol)
			
			append!(df_list,[df_temp])
		end
			
		spec_names = formatStrArr(collect(keys(sol_arr[1][1])))
		prepend!(spec_names,["timestamp"])
		df_dict = Dict()
		for i in 1:length(df_list)
			name = "Plot "*string(i)
			df_dict[Symbol(name)] =  ( collect(DataFrames.eachcol(df_list[i])), spec_names)
		end
	XLSX.writetable("sim_results.xlsx"; df_dict..., overwrite = true)	

	end
						md"""

 Download data:  $(DownloadButton(read("sim_results.xlsx"), "sim_results.xlsx")) 

"""
	else
		nothing
	end

		
end

# ╔═╡ afea37f1-70c2-4aae-94f6-34cf7c1d9f8e
begin
		sol2 = solve(ODEProblem(vf, cur_conc, (0.0,120.0),cur_rate, saveat=collect(0:120)));
	if importCheck == false
	
	df = DataFrame(sol2)
	nms = collect(keys(sol2(0)))
	prepend!(nms,[:timestamp])
	rename!(df,nms)
	nmsStr = formatStrArr(nms)
	filtered_names = []	
	filtered_names_inx = []
	
	if length(filter_list) > 0
		for i in 1:length(filter_list)
			name = filter_list[i]
				println(name)
				for j in 1:length(nmsStr)
					
					if occursin(name,nmsStr[j]) == true
						append!(filtered_names_inx,[j])
					end
				end
		end
	end
	sort!(unique!(filtered_names_inx))
	deleteat!(nmsStr,filtered_names_inx)
	
	finKeys = formatSymbArr(nmsStr)

	dfFin = select(df,finKeys)
		CSV.write("sim_res.csv", dfFin)
		
		
	# sol2.u
	# CSV.write("sim_res.csv", DataFrame(sol2), header = vcat([:timestamp], collect(keys(sol(0)))))
	else

	df = DataFrame(sol2)
	nms = collect(keys(sol2(0)))
	prepend!(nms,[:timestamp])
	rename!(df,nms)
		nmsStr = formatStrArr(nms)
	filtered_names = []	
	filtered_names_inx = []
	
	if length(filter_list) > 0
		for i in 1:length(filter_list)
			name = filter_list[i]
				println(name)
				for j in 1:length(nmsStr)
					
					if occursin(name,nmsStr[j]) == true
						append!(filtered_names_inx,[j])
					end
				end
		end
	end
	sort!(unique!(filtered_names_inx))
	deleteat!(nmsStr,filtered_names_inx)
	intersect!(nmsStr,graphKeys)
	
	finKeys = formatSymbArr(nmsStr)
	prepend!(finKeys,[:timestamp])
	
	dfFin = select(df,finKeys)
		CSV.write("sim_res.csv", dfFin)
	end

	md""" Download simulation data:  $(DownloadButton(read("sim_res.csv"), "sim_results.csv")) """
end

# ╔═╡ a141cd27-6ea0-4f73-80b5-72d8e5770ed4
begin
  if c != 0
		
		if filterInact == true || filterDeg == true
			timesteps = dfFin[!,"timestamp"]
			data = Matrix(dfFin[!,2:end])/1e3
			labels_new = names(dfFin)[2:end]
			plot(timesteps,data, label = reshape(labels_new, (1,length(labels_new))), linewidth = 3, xlabel = "Minutes",ylabel = "Solution Concentration (nM)")
			
		elseif combineCheck == false
		  tsteps = sol.t
		  # labels = [:G_deg]
			labels = isempty(graphKeySymb) ? snames(model) : graphKeySymb
		  plot(tsteps, [[sol(t)[l]/1e3 for t in tsteps] for l in labels], labels=hcat(String.(labels)...), linewidth=3, xlabel="Minutes", ylabel="Solution Concentration (nM)")

		else
			graphKeySymb2 = Symbol[]
			labels2 = isempty(graphKeySymb2) ? snames(model) : graphKeySymb2
			
			S_labels2 = formatStrArr(labels2)
			prepend!(S_labels2,["timestamp"])
			dfs2 = DataFrame(sol)
			rename!(dfs2,S_labels2)
			# names(dfs)[1][1]
			# dff2 = combineConc(dfs2)
			dff2 = combineConcNoDeg(dfs2)
			labels_new2 = names(dff2)[2:end]

			timesteps2 = dff2[!,"timestamp"]
			dff2[!,2:end]
			data2 = Matrix(dff2[!,2:end])/1e3
			plot(timesteps2,data2, label = reshape(labels_new2, (1,length(labels_new2))), linewidth = 3, xlabel = "Minutes",ylabel = "Solution Concentration (nM)")
			# Vector(labels_new)
		end

	end

end

# ╔═╡ 9625798a-67df-49e4-91ce-c7e23ed2a177
begin
	get_inds(d, names) = [k=>d[k] for k in names]
LabelledReactionNet{Number, Number}(model, get_inds(def_concs, snames(model)), get_inds(m_rates, tnames(model))) |> AffinityNet |> to_graphviz
end

# ╔═╡ e4100b5b-b255-48db-a989-016fa72f8da5
begin
sep_sym = "↦"
	
function split2(cat,sub,site::Int,on::T) where T
  catsym = first(cat)
  subsym = first(sub)
  in = Symbol(catsym,site,sep_sym,subsym)
  prod1 = Symbol(String(subsym)[1:site])
  prod2 = Symbol(String(subsym)[site+1:end])
  Open(LabelledReactionNet{T,Number}(unique((in=>0, cat,prod1=>0,prod2=>0)), ((Symbol(:split_,in),on),in=>(first(cat),prod1,prod2))));
end;

function split2(cat,sub,site::Int)
  in = Symbol(cat,site,sep_sym,sub)
  prod1 = Symbol(String(sub)[1:site])
  prod2 = Symbol(String(sub)[site+1:end])
  Open(LabelledPetriNet(unique((in, cat,prod1,prod2)), (Symbol(:split_,in),in=>(cat,prod1,prod2))));
end;

function bindunbind_multisite(in1, in2, site::Int, on::T, off::T) where T
  out = Symbol(first(in1),site, sep_sym,first(in2))
  Open(LabelledReactionNet{T,Number}(unique((in1, in2, out=>0)), ((Symbol(:bind,site,:_,first(in1), sep_sym,first(in2)),on),(first(in1),first(in2))=>out),
    ((Symbol(:unbind_,out),off),out=>(first(in1),first(in2)))))
end;

function bindunbind_multisite(in1, in2, site::Int)
  out = Symbol(in1,site, sep_sym,in2)
  Open(LabelledPetriNet(unique((in1, in2,out)), (Symbol(:bind,site,:_,in1, sep_sym,in2),(in1,in2)=>out),
                                                (Symbol(:unbind_,out),out=>(in1,in2))))
end;

enzXsubYZ = @relation (X, Xinact, Xdeg, YZ, Y, Z) where (X, Xinact, Xdeg, YZ, XYZ, Y, Z) begin
  bindXYZ(X, YZ, XYZ)
  splitXYZ(XYZ, X, Y, Z)
end

function enz_sub_split(rxns, cat1, sub, site)
	catsym = first(cat1)
  subsym = first(sub)
  
  catsub = Symbol(catsym, site, sep_sym, subsym)
  frag1sym = Symbol(String(subsym)[1:site])
  frag2sym = Symbol(String(subsym)[site+1:end])
    
  obtype = valtype(rates(apex(first(last(first(rxns))))))
  out = oapply(enzXsubYZ, Dict([:bindXYZ, :splitXYZ] .=> rxns[catsub]), Dict(
    :X=>ob(obtype, cat1),
    :Xinact=>ob(obtype, Symbol(catsym,:_inact)=>0),
    :Xdeg=>ob(obtype, Symbol(catsym,:_deg)=>0),
    :YZ=>ob(obtype, sub),
    :XYZ=>ob(obtype, catsub=>0),
    :Y=>ob(obtype, frag1sym=>0),
    :Z=>ob(obtype, frag2sym=>0)))
  println(sub," ",length(legs(out)))
	# bundle_legs(out, [[1,2,3], [4,5]])
	bundle_legs(out, [[1,2,3],[4],[5],[6]])
end;

function enz_sub_split(cat1::Symbol, sub::Symbol, site)
  catsym = cat1
  subsym = sub

  catsub = Symbol(catsym, site, sep_sym, subsym)
  frag1sym = Symbol(String(subsym)[1:site])
  frag2sym = Symbol(String(subsym)[site+1:end])

  out = oapply(enzXsubYZ, Dict(:bindXYZ=>bindunbind_multisite(catsym,subsym,site), :splitXYZ=>split2(catsym,subsym,site)), Dict(
    :X=>ob(cat1),
    :Xinact=>ob(Symbol(catsym,:_inact)),
    :Xdeg=>ob(Symbol(catsym,:_deg)),
    :YZ=>ob(subsym),
    :XYZ=>ob(catsub),
    :Y=>ob(frag1sym),
    :Z=>ob(frag2sym)))
  # bundle_legs(out, [[1,2,3], [4,5]])
  bundle_legs(out, [[1,2,3],[4],[5],[6]])
end;

function multisplit_generators(enzyme::Symbol, molecule::Symbol)
  gens = Dict{Symbol, Any}()
  # println((length(String(molecule))-1))
  for ii in 1:(length(String(molecule))-1)
    gens[Symbol(:cat, enzyme, ii, :sub, molecule)] = enz_sub_split(enzyme, molecule, ii)  
    frag1 = Symbol(String(molecule)[1:ii])
    frag2 = Symbol(String(molecule)[ii+1:end])
	# println(enzyme," ",frag1," ",frag2)
    if ii != 1
      merge!(gens,multisplit_generators(enzyme, frag1))
    end
    if ii != length(String(molecule))-1
      merge!(gens,multisplit_generators(enzyme, frag2))
    end
  end
  gens
end;

function multisplit_generators(enzymes::Array{Symbol}, molecule::Symbol)
  gens = Dict{Symbol, Any}()
  for enzyme1 in enzymes
    merge!(gens,multisplit_generators(enzyme1, molecule))
    for enzyme2 in enzymes
      if enzyme1==enzyme2
        gens[Symbol(:cat,enzyme1)] = enz(enzyme1)  
      else
        gens[Symbol(:cat,enzyme1,:cat,enzyme2)] = enz_enz(enzyme1,enzyme2)
      end
    end  
  end
  gens
end;
	
function gen_fragments(substrate::Symbol)
  frags = [substrate]
  for ii in 1:(length(String(substrate))-1)
    frag1 = Symbol(String(substrate)[1:ii])
    frag2 = Symbol(String(substrate)[ii+1:end])
    if ii != 1
      append!(frags,gen_fragments(frag1))
    end
    if ii != length(String(substrate))-1
      append!(frags,gen_fragments(frag2))
    end
  end
  frags
end;
	
function multisplit_uwd(enzymes::Array{Symbol}, substrate::Symbol) 
  rel = RelationDiagram{Symbol}(0)

  substrates = gen_fragments(substrate)
  chemicals = vcat(substrates, enzymes)

  subs = add_junctions!(rel, length(substrates), variable=substrates)
  enzs = add_junctions!(rel, length(enzymes), variable=enzymes)
  nsubs = length(subs)
  nenzs = length(enzs)

  catx = add_parts!(rel, :Box, nenzs, name=[Symbol("cat$i") for i in enzymes])
  add_parts!(rel, :Port, nenzs, junction=enzs, box=catx)

  for x in 1:nenzs
    for y in 1:nenzs
      if y != x
        catxy = add_part!(rel, :Box, name=Symbol("cat$(enzymes[x])cat$(enzymes[y])"))
        add_parts!(rel, :Port, 2, junction=[enzs[x], enzs[y]], box=catxy)
      end
    end
  end

  for x in 1:nenzs
    for y in 1:nsubs
      for n in 1:(length(String(substrates[y]))-1)
        catxy = add_part!(rel, :Box, name=Symbol("cat$(enzymes[x])$(n)sub$(substrates[y])"))
        add_parts!(rel, :Port, 2, junction=[enzs[x], subs[y]], box=catxy)
      end
    end
  end
  add_parts!(rel, :OuterPort, length(chemicals), outer_junction = vcat(subs, enzs))
  rel
end;
end

# ╔═╡ ef14041a-9753-456d-a242-d08dc0328507
begin
  # Initial Concentrations
  ABC = :ABC=>66000;
  AB = :AB=>0;
  BC = :BC=>0;
  A = :A=>0;
  B = :B=>0;
  subs2 = [:ABC]

  # Parameter Rates (units of pM and min)
  rxns2 = Dict(
    "K" => [inactivate(K, 7.494e-10)
           bindunbind(K, K, 7.814e-4, 3.867e-3)
           degrade(K, K, 2.265e-1)
           bindunbind(K, Kinact, 7.814e-4, 3.867e-3)
           degrade(K, Kinact, 2.265e-1)],
    "S" => [inactivate(S, 1.906e-2)
           bindunbind(S, S, 3.534e-7, 1.688e2)
           degrade(S, S, 1.433e1)
           bindunbind(S, Sinact, 3.534e-7, 1.688e2)
           degrade(S, Sinact, 1.433e1)],
    "L" => [inactivate(L, 7.810e-3)
           bindunbind(L, L, 1.000e-11, 7.440e3)
           degrade(L, L, 2.670e2)
           bindunbind(L, Linact, 1.000e-11, 7.440e3)
           degrade(L, Linact, 2.670e2)],
    "K1↦ABC" => [bindunbind_multisite(K, ABC, 1,9.668e-6,1.000e-2)
            split2(K, ABC, 1, 1.728e0)],
    "K2↦ABC" => [bindunbind_multisite(K, ABC, 2,9.668e-6,1.000e-2)
            split2(K, ABC, 2, 1.728e0)],
    "K1↦AB" => [bindunbind_multisite(K, AB, 1,9.668e-6,1.000e-2)
            split2(K, AB, 1, 1.728e0)],
    "K1↦BC" => [bindunbind_multisite(K, BC, 1,9.668e-6,1.000e-2)
            split2(K, BC, 1, 1.728e0)],
    "L1↦ABC" => [bindunbind_multisite(L, ABC, 1,9.668e-6,1.000e-2)
            split2(L, ABC, 1, 1.728e0)],
    "L2↦ABC" => [bindunbind_multisite(L, ABC, 2,9.668e-6,1.000e-2)
            split2(L, ABC, 2, 1.728e0)],
    "L1↦AB" => [bindunbind_multisite(L, AB, 1,9.668e-6,1.000e-2)
            split2(L, AB, 1, 1.728e0)],
    "L1↦BC" => [bindunbind_multisite(L, BC, 1,9.668e-6,1.000e-2)
            split2(L, BC, 1, 1.728e0)],
    "S1↦ABC" => [bindunbind_multisite(S, ABC, 1,9.668e-6,1.000e-2)
            split2(S, ABC, 1, 1.728e0)],
    "S2↦ABC" => [bindunbind_multisite(S, ABC, 2,9.668e-6,1.000e-2)
            split2(S, ABC, 2, 1.728e0)],
    "S1↦AB" => [bindunbind_multisite(S, AB, 1,9.668e-6,1.000e-2)
            split2(S, AB, 1, 1.728e0)],
    "S1↦BC" => [bindunbind_multisite(S, BC, 1,9.668e-6,1.000e-2)
            split2(S, BC, 1, 1.728e0)],
    "K↦S" => [bindunbind(K, S, 8.822e-4, 4.114e5)
            degrade(K, S, 9.000e-10)
            bindunbind(K, Sinact, 8.822e-4, 4.114e5)
            degrade(K, Sinact, 9.000e-10)],
    "K↦L" => [bindunbind(K, L, 1.756e-4, 3.729e4)
            degrade(K, L, 6.505e6)
            bindunbind(K, Linact, 1.756e-4, 3.729e4)
            degrade(K, Linact, 6.505e6)],
    "S↦K" => [bindunbind(S, K, 3.679e-4, 1.562e3)
            degrade(S, K, 4.410e2)
            bindunbind(S, Kinact, 3.679e-4, 1.562e3)
            degrade(S, Kinact, 4.410e2)],
    "S↦L" => [bindunbind(S, L, 1.000e-3, 5.000e2)
            degrade(S, L, 1.000e-7)
            bindunbind(S, Linact, 1.000e-3, 5.000e2)
            degrade(S, Linact, 1.000e-7)],
    "L↦K" => [bindunbind(L, K, 1.000e-3, 4.118e3)
            degrade(L, K, 3.234e1)
            bindunbind(L, Kinact, 1.000e-3, 4.118e3)
            degrade(L, Kinact, 3.234e1)],
    "L↦S" => [bindunbind(L, S, 1.056e-12, 5.000e2)
            degrade(L, S, 5.000e-1)
            bindunbind(L, Sinact, 1.056e-12, 5.000e2)
            degrade(L, Sinact, 5.000e-1)]
  );
  rxns2 = Dict(Symbol(k)=>v for (k,v) in rxns2)
  # define labels to reaction network mappings
  functor2(x) = oapply(x, Dict(
    :catK=>enz(rxns2, K),
    :catS=>enz(rxns2, S),
    :catL=>enz(rxns2, L),
    :catKcatS=>enz_enz(rxns2, K,S),
    :catKcatL=>enz_enz(rxns2, K,L),
    :catScatK=>enz_enz(rxns2, S,K),
    :catScatL=>enz_enz(rxns2, S,L),
    :catLcatK=>enz_enz(rxns2, L,K),
    :catLcatS=>enz_enz(rxns2, L,S),
    :catK1subABC=>enz_sub_split(rxns2, K,ABC,1),
    :catK2subABC=>enz_sub_split(rxns2, K,ABC,2),
    :catS1subABC=>enz_sub_split(rxns2, S,ABC,1),
    :catS2subABC=>enz_sub_split(rxns2, S,ABC,2),
    :catL1subABC=>enz_sub_split(rxns2, L,ABC,1),
    :catL2subABC=>enz_sub_split(rxns2, L,ABC,2),
    :catK1subAB=>enz_sub_split(rxns2, K,AB,1),
    :catS1subAB=>enz_sub_split(rxns2, S,AB,1),
    :catL1subAB=>enz_sub_split(rxns2, L,AB,1),
    :catK1subBC=>enz_sub_split(rxns2, K,BC,1),
	:catS1subBC=>enz_sub_split(rxns2, S,BC,1),
	:catL1subBC=>enz_sub_split(rxns2, L,BC,1)
	));
  lfunctor2(x) = oapply(x, multisplit_generators([:K,:L],:ABC))
  def_model2 = apex(functor2(multisplit_uwd([:K,:L],:ABC)))
  def_rates2 = rates(def_model2)
  def_concs2 = Dict(c=>concentrations(def_model2)[c] for c in snames(def_model2))
  def_concs2[:Spike] = Spike[2]
  def_concs2[:Spike_deg] = 0
  def_concs2[Symbol("K↦Spike")] = 0
  def_concs2[Symbol("S↦Spike")] = 0
  def_concs2[Symbol("L↦Spike")] = 0
  nothing
end

# ╔═╡ ecdc5f61-6041-42ef-819c-1d83c062c8e3
Graph(def_model2)

# ╔═╡ 7f6b0fe8-8ae1-4a42-87ea-2d5d3ae95181
multisplit_uwd([:K,:L], :ABC) |> lfunctor2 |> apex |> Graph

# ╔═╡ 41fa1014-0d0a-4b30-9452-c5a1fc8e58b5
md""" Upload a rate data file to use:  $(@bind user_csv2 FilePicker()) """

# ╔═╡ 553bce96-a29f-4a77-a193-8005914f4bfa
md"""Use default rates (override uploaded rates)? $(@bind defbox2 CheckBox(false))"""

# ╔═╡ e5d3b132-baca-4cdb-aeb0-09272610ed6f
begin
	local k_rates, status
	try
		defbox2 && error("default rate override")
		k_rates = UInt8.(user_csv2["data"]) |> IOBuffer |> JSON.parse
		status = md"""Using rates from $(user_csv2["name"])"""
	catch e
		k_rates = def_rates2
		status = md"""No file selected, reverting to default rates"""
	end
	m_rates2 = Dict(Symbol(k)=>k_rates[k] for k in keys(k_rates))
	valid_enz_sub2 = Set{Symbol}()
	for e in cats
		for s in cats
      if Symbol("bind_$e↦$s") ∈ keys(m_rates2)
				push!(valid_enz_sub2, e)
				push!(valid_enz_sub2, s)
			end
		end
		for s in gen_fragments(subs2[1])
			for i in 1:(length(String(s))-1)
      if Symbol("bind$(i)_$(e)↦$(s)") ∈ keys(m_rates2)
				push!(valid_enz_sub2, e)
				push!(valid_enz_sub2, s)
			end
			end
		end
	end
	status
end

# ╔═╡ 240a8494-158f-4db2-9d0d-c141b50dcd5d
md"""
#### Enzyme Diagram
Enzymes

"""

# ╔═╡ 48921b9a-d3fd-4d03-9170-a414797d8dee
begin
#TODO
#Find a more compact way to represent these checkboxes while still dynamically #determining which to display

	local boxes = Dict(:K=>md"K $(@bind kbox2 CheckBox(:K ∈ valid_enz_sub2))",
				 :S=>md"S $(@bind sbox2 CheckBox(:S ∈ valid_enz_sub2))",
				 :L=>md"L $(@bind Lbox2 CheckBox(:L ∈ valid_enz_sub2 && false))")
md"""
$(:K ∈ valid_enz_sub2 ? boxes[:K] : md"")
$(:S ∈ valid_enz_sub2 ? boxes[:S] : md"")
$(:L ∈ valid_enz_sub2 ? boxes[:L] : md"")
"""
end

# ╔═╡ 693397a9-3c80-41e3-b618-fa9d5ecb42b7
md"""Substrates"""

# ╔═╡ d5a956d0-b443-4ee2-9045-14146012a435
begin

	local boxes = Dict(:G=>md"G $(@bind gbox2 CheckBox(:G ∈ valid_enz_sub2))",
				 :E=>md"E $(@bind ebox2 CheckBox(:E ∈ valid_enz_sub2 && false))",
				 :ABC=>md"ABC $(@bind abcbox2 CheckBox(:ABC ∈ valid_enz_sub2))",
				 :AB=>md"AB $(@bind abbox2 CheckBox(:AB ∈ valid_enz_sub2))",
				 :BC=>md"BC $(@bind bcbox2 CheckBox(:BC ∈ valid_enz_sub2))")
md"""
$(:G ∈ valid_enz_sub2 ? boxes[:G] : md"")
$(:E ∈ valid_enz_sub2 ? boxes[:E] : md"")
$(:ABC ∈ valid_enz_sub2 ? boxes[:ABC] : md"")
$(:AB ∈ valid_enz_sub2 ? boxes[:AB] : md"")
$(:BC ∈ valid_enz_sub2 ? boxes[:BC] : md"")
"""
end

# ╔═╡ 2012c352-3bff-4c84-bc2e-d83b22d95349
begin
inp2 = [kbox2, sbox2, Lbox2];
inpSymb2 = [:K, :S, :L];
inpF2 = Symbol[];
	for i = 1:length(inp2)
		if inp2[i]
			push!(inpF2,inpSymb2[i]);
		end

	end

outp2 = [gbox2, ebox2, abcbox2, abbox2, bcbox2];
outSymb2 = [:G, :E, :ABC, :AB, :BC];
outF2 = Symbol[];
	for i = 1:length(outp2)
		if outp2[i]
			push!(outF2,outSymb2[i]);
		end

	end

end


# ╔═╡ 6b819b32-3601-48cc-9ff4-1df3e88e8034
uwd2 = multisplit_uwd(inpF2, :ABC);

# ╔═╡ 0858af84-9c5c-43b4-b5c9-fec1f1c4dba4
begin
  display_uwd(uwd2)
end

# ╔═╡ ab8a18b1-4975-4f4a-994b-b08773489baf
begin 
model2 = uwd2 |> lfunctor2 |> apex;
  r2_1 = join(["'$k': $(m_rates2[k])" for k in tnames(model2)], ", ");

  r2_2 = join(["'$k': $(def_concs2[k])" for k in snames(model2)], ", ");
  #== Basic Slider Template
<li><input
  class='slider'
  type='range'
  min='-10'
  step='0.01'
  max='10'
  value='0'
  oninput='this.nextElementSibling.value=(10**this.value).toExponential(2)'
  ></input>
<output>1</output></li>
==#
c2 = vcat(["$(m_rates2[k])" for k in tnames(model2)], ["$(def_concs2[k])" for k in snames(model2)]);
	
form_vals2 = HTML("""
<bond def = "c2">
<form>
  <div id ="myDIV" style='height:500px;overflow:scroll'>
  <h4>Adjust Rates</h4>
  <table>
  </table>
  </div>
 <div id ="concDiv" style='height:500px;overflow:scroll'>
		<h4>Adjust Concentrations</h4>
		<table id = "table2">
		</table>
	</div>
 <div id ="bindingAff" style='height:500px;overflow:scroll'>
		<h4>Binding Affinities</h4>
		<table id = "bindingAffTable">
		</table>
	</div>
</form>
</bond>
<style>
#myDIV {
  width: 100%;
  padding: 50px 0;
  text-align: center;
 background-color: #F9F6E5;
  margin-top: 20px;
}
#concDiv {
  width: 100%;
  padding: 50px 0;
  text-align: center;
  background-color:  #F9F6E5;
  margin-top: 20px;
  display: none
}
#bindingAff {
  width: 100%;
  padding: 50px 0;
  text-align: center;
  background-color:  #F9F6E5;
  margin-top: 20px;
  display: none
}
.button {
  background-color: #003057;
  border: none;
  color: white;
  padding: 15px 32px;
  text-align: center;
  text-decoration: none;
  display: inline-block;
  font-size: 16px;
  margin: 4px 2px;
  cursor: pointer;
}
.button:focus {
    background-color:#BFB37C;
}
.slider {
  -webkit-appearance: none;
  width: 100%;
  height: 7px;
  border-radius: 5px;
  background: #d3d3d3;
  outline: none;
  opacity: 0.7;
  -webkit-transition: .2s;
  transition: opacity .2s;
}
.slider::-webkit-slider-thumb {
  -webkit-appearance: none;
  appearance: none;
  width: 20px;
  height: 20px;
  border-radius: 50%;
  background:  #BFB37C;
  cursor: pointer;
}
.slider::-moz-range-thumb {
  width: 20px;
  height: 20px;
  border-radius: 50%;
  background: #003057;
  cursor: pointer;
}
</style>
<script>
//`currentScript` is the current script tag - we use it to select elements//
function hideShowRate() {
	var x = document.getElementById("myDIV");
	var y = document.getElementById("concDiv");
	var z = document.getElementById("bindingAff");
	z.style.display = "none";
	x.style.display = "block";
	y.style.display = "none"
}
function hideShowConc() {
	var x = document.getElementById("concDiv");
	var y = document.getElementById("myDIV");
	var z = document.getElementById("bindingAff");
	x.style.display = "block";
	y.style.display = "none"
	z.style.display = "none"
}
function hideShowBind() {
	var x = document.getElementById("bindingAff");
	var y = document.getElementById("myDIV");
	var z = document.getElementById("concDiv");
	x.style.display = "block";
	y.style.display = "none"
	z.style.display = "none"
	generateBindingTable()
}
function returnStrings(s,arr){
		let finalArr = []
		for (let item in arr){
			if(arr[item].includes(s) == true){
				finalArr.push(arr[item])
				}
			}
		return finalArr
		}
function generateBindingTable(){
		let keyNames = Object.keys(rates)
		let objVals = Object.values(rates)
		let bindArr = []
		let unbindArr = []
		for (let i in keyNames){
		   let nm = keyNames[i]
			if(nm.substring(0,4) == 'bind'){
				bindArr.push(nm)
		}else if(nm.substring(0,6) == 'unbind'){
				unbindArr.push(nm)
			}
			}
		let kdArr = []; //array of binding affinities
		let kdNames = []; // array of binding affinity names
			//grab values from sliders
	    let x = form.getElementsByClassName('sliderval');
		let xarr = Array.from(x, (v,_)=>{return v.value})
		//Loop over bind and unbind arrays to do calculations
		for(let i in keyNames){
			if(keyNames[i].substring(0,4) == 'bind'){
			let bind = keyNames[i]
			for(let j in keyNames){
		    if(keyNames[j].substring(0,6) == 'unbind'){
				let unbind = keyNames[j]
					if(bind.substring(5,bind.length) == unbind.substring(7,unbind.length)){
								kdArr.push(xarr[j]/xarr[i])
						kdNames.push('Kd_'+bind.substring(5,bind.length))
		//Create name of deg
		let nm = bind.substring(4,bind.length)
		let deg_name = 'deg' + nm
			for(let k in keyNames){
				if(keyNames[k]== deg_name){
						let kcat = keyNames[k]
						let unbindval = parseFloat(xarr[j]);
						let bindval = parseFloat(xarr[i])
						let catval = parseFloat(xarr[k])
						let k_m = (catval+unbindval)/bindval;
						kdArr.push(k_m)
						kdNames.push('Km_'+bind.substring(5,bind.length))
						kdArr.push(catval/k_m)
						kdNames.push('Cat_efficieny_'+bind.substring(5,bind.length))
							}
						}}}
					}
				}
			}
		const list3 = form.querySelector('#bindingAffTable')
		while(list3.rows.length > 0) {
				  list3.deleteRow(0);
				}
		var header = document.createElement('thead')
 		var toprow = document.createElement('tr');
		var toprowItem = document.createElement('th');
		var toprowItem2 = document.createElement('th');
		var toprowItem3 = document.createElement('th');
		var toprowItem4 = document.createElement('th');
		toprowItem.innerText = 'Name';
		toprowItem2.innerText = 'Kd';
		toprowItem3.innerText = 'Km';
		toprowItem4.innerText = 'Catalytic Efficiency';
		header.appendChild(toprow)
		//toprow.appendChild(header);
		toprow.appendChild(toprowItem);
		toprow.appendChild(toprowItem2);
		toprow.appendChild(toprowItem3);
		toprow.appendChild(toprowItem4);
		list3.appendChild(toprow);
		for (let i = 0; i < kdNames.length/3; i++){
		  var item = document.createElement('tr');
		  var label = document.createElement('td');
		  label.innerText = kdNames[i*3].substring(3)
		  var label2 = document.createElement('td');
		  var label3 = document.createElement('td');
		  var label4 = document.createElement('td');
		  label2.innerText =  kdArr[i*3].toExponential(2)
		label3.innerText = kdArr[i*3+1].toExponential(2) // 'placeholder'
		label4.innerText = kdArr[i*3+2].toExponential(2) // 'placeholder'
		  item.appendChild(header)
		  item.appendChild(label)
          item.appendChild(label2)
		  item.appendChild(label3)
		  item.appendChild(label4)
		  list3.appendChild(item)
		}
		}
const form = currentScript.parentElement.querySelector('form')
const list = form.querySelector('table')
var rates = {$r2_1};
var conc = {$r2_2};
for ( var r in rates ){
  console.log(r);
  var item = document.createElement('tr');
  var label = document.createElement('th');
  label.innerText = r
  var slider_box = document.createElement('td');
  var slider_val_box = document.createElement('td');
  var slider = document.createElement('input');
  var slider_val = document.createElement('input');
  slider_val.setAttribute('type', 'text');
  slider_val.setAttribute('class', 'sliderval');
  slider_val.value = rates[r].toExponential(2);
  slider.setAttribute('type', 'range');
  slider.setAttribute('class', 'slider');
  slider.setAttribute('min', '-9.99');
  slider.setAttribute('max', '10');
  slider.setAttribute('value', '0.0');
  slider.setAttribute('step', '0.01');
  slider.setAttribute('oninput', `this.parentElement.nextElementSibling.children[0].value=((10**this.value)*\${rates[r]}).toExponential(2)`);
  slider_box.appendChild(slider)
  slider_val_box.appendChild(slider_val)
  item.appendChild(label)
  item.appendChild(slider_box)
  item.appendChild(slider_val_box)
  list.appendChild(item)
}
const list2 = form.querySelector('#table2')
for ( var r in conc ){
  console.log(r);
  var item = document.createElement('tr');
  var label = document.createElement('th');
  label.innerText = r
  var slider_box = document.createElement('td');
  var slider_val_box = document.createElement('td');
  var slider = document.createElement('input');
  var slider_val = document.createElement('input');
  slider_val.setAttribute('type', 'text');
  slider_val.setAttribute('class', 'sliderval');
  slider_val.value = conc[r].toExponential(2);
  slider.setAttribute('type', 'range');
  slider.setAttribute('class', 'slider');
  slider.setAttribute('min', '-9.99');
  slider.setAttribute('max', '10');
  slider.setAttribute('value', '0.0');
  slider.setAttribute('step', '0.01');
  console.log('concentration')
  console.log(conc[r]);
  let test = conc[r];
  if(test == 0){
		test = .000001
		}
  slider.setAttribute('oninput', `this.parentElement.nextElementSibling.children[0].value=((10**this.value)*\${test}).toExponential(2)`);
  slider_box.appendChild(slider)
  slider_val_box.appendChild(slider_val)
  item.appendChild(label)
  item.appendChild(slider_box)
  item.appendChild(slider_val_box)
  list2.appendChild(item)
}
var x = form.getElementsByClassName('sliderval');
function onsubmit(){
  // We send the value back to Julia //
  form.value = Array.from(x, (v,_)=>{return v.value})
  form.dispatchEvent(new CustomEvent('input'))
  console.log(form.value)
}
var b = document.createElement('input');
b.setAttribute('type', 'button');
b.setAttribute('class', 'button');
b.value = 'Update Plot';
b.addEventListener('click', function() {onsubmit();console.log('hello from button')})
form.appendChild(b)
onsubmit()
let sp2 = document.getElementById("myDIV")
var d = document.createElement('input')
d.setAttribute('type','button')
d.setAttribute('class','button')
d.value = 'Rates'
d.addEventListener('click', function() {hideShowRate();})
form.insertBefore(d,sp2)
var e = document.createElement('input')
e.setAttribute('type','button')
e.setAttribute('class','button')
e.value = 'Concentrations'
e.addEventListener('click', function() {hideShowConc();})
form.insertBefore(e,sp2)
var f = document.createElement('input')
f.setAttribute('type','button')
f.setAttribute('class','button')
f.value = 'Binding Affinities'
f.addEventListener('click', function() {hideShowBind();})
form.insertBefore(f,sp2)
</script>
""");
end

# ╔═╡ 79064f4a-f4ad-4837-898a-d7a7e03f04da
begin
cur_rate2 = Dict(tnames(model2)[i]=>parse(Float64, c2[i]) for i in 1:length(tnames(model2)))
  # cur_conc = concentrations(model)
  cur_conc_prep2 =  Dict(snames(model2)[i]=>parse(Float64, c2[i+length(tnames(model2))]) for i in 1:length(snames(model2)))

	keysa2 = keys(cur_conc_prep2)
	keyArr2_1 = Symbol[]
	keyArrStr2 = String[]
	valArr2 = Float64[]
	for key in keys(cur_conc_prep2)
		push!(keyArr2_1,key)
		push!(valArr2,cur_conc_prep2[key])
		push!(keyArrStr2,String(key))
	end

  cur_conc2 = @LArray valArr2 Tuple(keyArr2_1)
  vf2 = vectorfield(model2);
sol2_1 = solve(ODEProblem(vf2, cur_conc2, (0.0,120.0),cur_rate2));
  nothing
end

# ╔═╡ 522d5f7b-7774-40f3-aa66-6dba2b48783a
length(c2)

# ╔═╡ 31a3d715-4a80-4696-96c3-8e5152797f9f
md"""### Iterate over multiple species concentrations $(@bind multCheck2 CheckBox(false))
#### Select species to iterate: 
$(@bind veg2 Select(keyArrStr))

Starting Concentration: $(@bind initConc2 TextField()) \
End Concentration: $(@bind endConc2 TextField()) \
number of steps: $(@bind stepsConc2 NumberField(1:10, default=5)) \

Choose spacing type: $(@bind spacing2 Select(["Logarithmic","Linear"])) \
Display combined concentrations excluding degraded? $(@bind combineCheckDeg2 CheckBox(false)) \
Use constant scale for axis? $(@bind scaleCheck2 CheckBox(false))

"""

# ╔═╡ 7a4fbc79-9221-472f-acfa-d32249c8e828
begin
if multCheck2
	if spacing2 == "Linear"
		conList2 = range(parse(Float64,initConc2), parse(Float64,endConc2), length = stepsConc2)
		
	elseif spacing2 == "Logarithmic"
		conList2 = logrange(initConc2,endConc2,stepsConc2);
		# println(conList)
	else
		conList2 = ["Empty"]
	end
	multiConcArr2 = createMultiArr(cur_conc2,conList2,veg2)
	sol_arr2 = []
	plotArr2 = []
	plotNumArr2 = []
	
	for item in multiConcArr2
		tempsol = solve(ODEProblem(vf2, item, (0.0,120.0),cur_rate2), saveat = 1.0);
		append!(sol_arr2,[tempsol])	
	end
		
	numPlots2 = length(multiConcArr2)
	for p in 1:length(sol_arr2)
		tsol = sol_arr2[p]
		concName = conList2[p]
		if combineCheckDeg2 == false
			append!(plotArr2,[createPlot(tsol,model2,concName,p)])
		else
			append!(plotArr2,[createPlotDF(tsol,model2,concName,p)])	
		end
		
	end
	
	# labels_mult = isempty(graphKeySymb) ? snames(model) : graphKeySymb
	labels_mult2 = collect(keys(multiConcArr2[1]))
	for i in 1:length(conList2)
		append!(plotNumArr2,[string("Plot ",string(i))])
	end
		if scaleCheck2 == false

			plot(plotArr2..., size = (600, 400*length(multiConcArr2)),layout = (length(multiConcArr2),1), legend = false, bottom_margin = 10mm)
		else
			maxlim2 = 1.0
			for j in 1:length(sol_arr2)
				
				if maximum(sol_arr2[j])>maxlim2
					maxlim2 = maximum(sol_arr2[j])
				end
			end
				plot(plotArr2..., size = (600, 400*length(multiConcArr2)),layout = (length(multiConcArr2),1), legend = false, bottom_margin = 10mm, ylims = (0,maxlim2/1e3))	
		end
	else
		plotNumArr2 = ["Empty"];
		nothing;
		
end
	
end

# ╔═╡ 59ed7708-3c09-431e-8a09-1e94b0b90e2a
begin
if multCheck2
keyArr3_2 = ["A","B"];
graphKeys3_2 = []
graphKeyVals3_2 = HTML("""
<bond def = "graphKeys3_2">
<form>

  <div id ="myDIV3" style='height:300px;overflow:scroll'>

  <h4>Select Species to plot</h4>
			<p>This plot will display how selected species change over iterations. Limit of 6 species. Iteration is indicated by color and species is indicated by marker shape</p>
			
  <select  id="graphConc" multiple = "multiple" size = "10">

  </select>



  </div>



</form>
</bond>
<style>

#myDIV2 {
  width: 100%;
  padding: 10px 0;
  text-align: center;
  background-color: #F9F6E5;
  margin-top: 10px;
}


</style>



<script>
//`currentScript` is the current script tag - we use it to select elements//


const form = currentScript.parentElement.querySelector('form')
const list =  form.querySelector('select')
console.log(list)
var concList = $keyArrStr2;
console.log(concList)

for ( var key in concList ){
  console.log(concList[key]);
  //var item = document.createElement('tr');
  //var label = document.createElement('th');
  var item = document.createElement('option');
  item.value = concList[key];
  item.label = concList[key];
  list.appendChild(item)

}





var x = form.getElementsByClassName('selector');
function onsubmit(){
		console.log('onsubmit called')
  // We send the value back to Julia //
		console.log(x)


    var selected = [];
    for (var option of document.getElementById('graphConc').options)
    {
        if (option.selected) {
            selected.push(option.value);
        }
    }
 form.value = selected
  form.dispatchEvent(new CustomEvent('input'))

  console.log(form.value)
    console.log(selected)
}
var b = document.createElement('input');
b.setAttribute('type', 'button');
b.setAttribute('class', 'button');
b.value = 'Plot species';
b.addEventListener('click', function() {onsubmit();console.log('hello from button')})
form.appendChild(b)
onsubmit()




</script>
""");
	else
		nothing;
	end
end

# ╔═╡ cf516855-40fa-44d0-8d63-c145ef4577d0
begin 
	df_array2 = [];
if multCheck2
		
		#Convert solutions to dataframe
		str_labels_mult2 = formatStrArr(labels_mult2)
		numSpecs2 = length(graphKeys3_2)
		nplots2 = length(sol_arr2)
		prepend!(str_labels_mult2,["timestamp"])
		timestamps2 = []
	for i in 1:length(sol_arr2)
			
		tempdf2 = DataFrame(sol_arr2[i])
		rename!(tempdf2,str_labels_mult2)
		timestamps2 = tempdf2[!, "timestamp"]
		append!(df_array2,[select(tempdf2,graphKeys3_2)])
	end
	df_plot2 = df_array2[1]
	# timestamps = sol_arr[1][!,:timestamp]
	for i in 2:length(sol_arr2)
		df_plot2 = hcat(df_plot2,df_array2[i], makeunique=true)
	end
		# println(df_plot)
		
		# mat_data = Matrix(df_plot)
		# plotnames = names(df_plot)
		# plotlabel = reshape(plotnames, (1,length(plotnames)))
		# plot(timestamps,mat_data, label = plotLabel, linewidth = 3, xlabel = "Minutes",ylabel = "Solution Concentration (nM)")
	else
		nothing;
	end
end

# ╔═╡ 6734cf5a-a3dc-4de3-9f33-45e149a27ebc
begin 
	if multCheck2 == true && length(df_array2)>0 
generatePlotDF2(df_plot2,timestamps2,numSpecs2,nplots2)
	else
		nothing;
	end
end

# ╔═╡ d9a98759-e596-4c81-98f5-2d4efdc964f6
md"""#### Choose export mode $(@bind exportMode2 Select(["Export all plots", "Export specific plot"]))
##### Select plot for export: 
$(@bind conListSel2 Select(formatStrArr(plotNumArr2)))
"""

# ╔═╡ 20c7902d-ebbf-49d1-a103-fb7524489524
begin
	if multCheck2
	if exportMode2 == "Export specific plot"
		ind2 = parse(Int8,conListSel2[end])
		println("EXporting specific Plot")
		sol_sel2 = sol_arr2[ind2][1]
		df_sel2 = DataFrame(sol_arr2[ind2])
		nms_sel2 = collect(keys(sol_sel2))
		prepend!(nms_sel2,[:timestamp])
		rename!(df_sel2,nms_sel2)
		nmsStr_sel2 = formatStrArr(nms_sel2)
		filtered_names_sel2 = []	
		filtered_names_inx_sel2 = []

		if length(filter_list) > 0
			for i in 1:length(filter_list)
				name = filter_list[i]
					
					for j in 1:length(nmsStr_sel2)

						if occursin(name,nmsStr_sel2[j]) == true
							append!(filtered_names_inx_sel2,[j])
						end
					end
			end
		end
		sort!(unique!(filtered_names_inx_sel2))
		deleteat!(nmsStr_sel2,filtered_names_inx_sel2)

		finKeys_sel2 = formatSymbArr(nmsStr_sel2)

		dfFin_sel2 = select(df_sel2,finKeys_sel2)
		# CSV.write("sim_res_sel.csv", dfFin_sel)
		XLSX.writetable("sim_results2.xlsx", collect(DataFrames.eachcol(dfFin_sel2)), DataFrames.names(dfFin_sel2);overwrite = true)

		
	elseif exportMode2 == "Export all plots"
		println("Export all plots")
		# Use the xlsx package to create an excel sheet with multiple tabs. Then can export dataframes
		#First convert solution into dataframes
		df_list2 = []

		for sltn in 1:length(plotNumArr2)
			println(sltn)
			temp_sol = sol_arr2[sltn]
			# println(temp_sol)
			df_temp = DataFrame(temp_sol)
			
			append!(df_list2,[df_temp])
		end
			
		spec_names2 = formatStrArr(collect(keys(sol_arr2[1][1])))
		prepend!(spec_names2,["timestamp"])
		df_dict2 = Dict()
		for i in 1:length(df_list2)
			name = "Plot "*string(i)
			df_dict2[Symbol(name)] =  ( collect(DataFrames.eachcol(df_list2[i])), spec_names2)
		end
	XLSX.writetable("sim_results2.xlsx"; df_dict2..., overwrite = true)	

	end
						md"""

 Download data:  $(DownloadButton(read("sim_results2.xlsx"), "sim_results2.xlsx")) 

"""
	else
		nothing
	end

		
end

# ╔═╡ 8c14d4a9-e2e7-4b18-bbba-a290872c2ca6
md"""Display total concentrations? $(@bind combineCheck2 CheckBox(false)) 

"""

# ╔═╡ bf7ea880-5dff-419d-ba69-b4271b04bbe2
begin

keyArr2_2 = ["A","B"];
graphKeys2 = []
graphKeyVals2 = HTML("""
<bond def = "graphKeys2">
<form>

  <div id ="myDIV2" style='height:300px;overflow:scroll'>

  <h4>Select Variables to Export</h4>
  <select  id="graphConc" multiple = "multiple" size = "10">

  </select>



  </div>



</form>
</bond>
<style>

#myDIV2 {
  width: 100%;
  padding: 10px 0;
  text-align: center;
  background-color: #F9F6E5;
  margin-top: 10px;
}


</style>



<script>
//`currentScript` is the current script tag - we use it to select elements//


const form = currentScript.parentElement.querySelector('form')
const list =  form.querySelector('select')
console.log(list)
var concList = $keyArrStr2;
console.log(concList)

for ( var key in concList ){
  console.log(concList[key]);
  //var item = document.createElement('tr');
  //var label = document.createElement('th');
  var item = document.createElement('option');
  item.value = concList[key];
  item.label = concList[key];
  list.appendChild(item)

}





var x = form.getElementsByClassName('selector');
function onsubmit(){
		console.log('onsubmit called')
  // We send the value back to Julia //
		console.log(x)


    var selected = [];
    for (var option of document.getElementById('graphConc').options)
    {
        if (option.selected) {
            selected.push(option.value);
        }
    }
 form.value = selected
  form.dispatchEvent(new CustomEvent('input'))

  console.log(form.value)
    console.log(selected)
}
var b = document.createElement('input');
b.setAttribute('type', 'button');
b.setAttribute('class', 'button');
b.value = 'Update variable list';
b.addEventListener('click', function() {onsubmit();console.log('hello from button')})
form.appendChild(b)
onsubmit()




</script>
""");

end

# ╔═╡ 5dc80a12-9265-4e5d-92b3-be8883f65ee7
md"""Filter Inact? $(@bind filterInact2 CheckBox(false)) \
Filter degraded? $(@bind filterDeg2 CheckBox(false))
"""


# ╔═╡ def4bd87-532b-4552-8229-620a27c04641
begin
	filter_list2 = [];
	
if filterInact2 == true || filterDeg2 == true
	if filterInact2 == true
			append!(filter_list2,["inact"]);
	end
	if filterDeg2 == true
			append!(filter_list2,["deg"]);
	end
		
end
	nothing;
end

# ╔═╡ 1b9fe2e0-6a17-4f6b-bb23-f371e261e86e
md"""Export only selected variables? $(@bind importCheck2 CheckBox(false))"""

# ╔═╡ 3eadc95c-b063-4ca6-ab5f-20acf5c56c2e
begin
		sol2_2 = solve(ODEProblem(vf2, cur_conc2, (0.0,120.0),cur_rate2, saveat=collect(0:120)));
	if importCheck == false
	
	df2 = DataFrame(sol2_2)
	nms2 = collect(keys(sol2_2(0)))
	prepend!(nms2,[:timestamp])
	rename!(df2,nms2)
	nmsStr2 = formatStrArr(nms2)
	filtered_names2 = []	
	filtered_names_inx2 = []
	
	if length(filter_list2) > 0
		for i in 1:length(filter_list2)
			name = filter_list2[i]
				println(name)
				for j in 1:length(nmsStr2)
					
					if occursin(name,nmsStr2[j]) == true
						append!(filtered_names_inx2,[j])
					end
				end
		end
	end
	sort!(unique!(filtered_names_inx2))
	deleteat!(nmsStr2,filtered_names_inx2)
	
	finKeys2 = formatSymbArr(nmsStr2)

	dfFin2 = select(df2,finKeys2)
		CSV.write("sim_res2.csv", dfFin2)
		
		
	# sol2.u
	# CSV.write("sim_res.csv", DataFrame(sol2), header = vcat([:timestamp], collect(keys(sol(0)))))
	else

	df2 = DataFrame(sol2_2)
	nms2 = collect(keys(sol2_2(0)))
	prepend!(nms2,[:timestamp])
	rename!(df2,nms2)
		nmsStr2 = formatStrArr(nms2)
	filtered_names2 = []	
	filtered_names_inx2 = []
	
	if length(filter_list) > 0
		for i in 1:length(filter_list)
			name = filter_list[i]
				println(name)
				for j in 1:length(nmsStr2)
					
					if occursin(name,nmsStr2[j]) == true
						append!(filtered_names_inx2,[j])
					end
				end
		end
	end
	sort!(unique!(filtered_names_inx2))
	deleteat!(nmsStr2,filtered_names_inx2)
	intersect!(nmsStr2,graphKeys2)
	
	finKeys2 = formatSymbArr(nmsStr2)
	prepend!(finKeys2,[:timestamp])
	
	dfFin2 = select(df2,finKeys2)
		CSV.write("sim_res2.csv", dfFin2)
	end

	md""" Download simulation data:  $(DownloadButton(read("sim_res2.csv"), "sim_results2.csv")) """
end

# ╔═╡ 9db5401c-8501-45b5-8d8e-070afcdfd867
begin
  if c2 != 0
		
		if filterInact2 == true || filterDeg2 == true
			timesteps2b = dfFin2[!,"timestamp"]
			data2b = Matrix(dfFin2[!,2:end])/1e3
			labels_new2b = names(dfFin2)[2:end]
			plot(timesteps2b,data2b, label = reshape(labels_new2b, (1,length(labels_new2b))), linewidth = 3, xlabel = "Minutes",ylabel = "Solution Concentration (nM)")
			
		elseif combineCheckDeg2 == false
		  tsteps2 = sol2_1.t
		  # labels = [:G_deg]
			labels2b = isempty(graphKeySymb) ? snames(model2) : graphKeySymb
		  plot(tsteps2, [[sol2_1(t)[l]/1e3 for t in tsteps2] for l in labels2b], labels=hcat(String.(labels2b)...), linewidth=3, xlabel="Minutes", ylabel="Solution Concentration (nM)")

		else
			graphKeySymb2b = Symbol[]
			labels2b = isempty(graphKeySymb2b) ? snames(model2) : graphKeySymb2b
			
			S_labels2_2 = formatStrArr(labels2b)
			prepend!(S_labels2_2,["timestamp"])
			dfs2_2 = DataFrame(sol2_1)
			rename!(dfs2_2,S_labels2_2)
			# names(dfs)[1][1]
			# dff2 = combineConc(dfs2)
			dff2_2 = combineConcNoDeg(dfs2_2)
			labels_new2_2 = names(dff2_2)[2:end]

			timesteps2_2 = dff2_2[!,"timestamp"]
			dff2_2[!,2:end]
			data2_2 = Matrix(dff2_2[!,2:end])/1e3
			plot(timesteps2_2,data2_2, label = reshape(labels_new2_2, (1,length(labels_new2_2))), linewidth = 3, xlabel = "Minutes",ylabel = "Solution Concentration (nM)")
			# Vector(labels_new)
		end

	end

end

# ╔═╡ 4886158c-8dca-4d70-9903-7287624983ef
begin

graphKeySymb2_2 = Symbol[]
for item in graphKeys2
		push!(graphKeySymb2_2,Symbol(item))

end
end

# ╔═╡ 33a0fc85-1069-4355-b2a6-6165f91d24e7
begin
	get_inds2(d, names) = [k=>d[k] for k in names]
LabelledReactionNet{Number, Number}(model2, get_inds2(def_concs2, snames(model2)), get_inds2(m_rates2, tnames(model2))) |> AffinityNet |> to_graphviz
end

# ╔═╡ ca18c7c8-6c58-4bd5-a92b-ed044c012603
valid_enz_sub2

# ╔═╡ bced5a4c-9e43-4910-b233-738675a0f887
gen_fragments(subs2[1])

# ╔═╡ Cell order:
# ╠═4c9c24cc-b865-4825-a841-f717120d27d2
# ╠═db1cf041-8f29-4433-a8b7-da5ea4828d52
# ╠═32c8703f-6aa3-46be-a91b-ff36225d6bd8
# ╠═aa5b099e-771f-42ff-a1b1-d22f0e615395
# ╠═178e764e-e239-4689-bb2f-4993b7755724
# ╠═563cf0a2-80e0-4bc2-8f6f-6a47cb2112af
# ╠═2d89b8e5-31a0-402c-b95c-87494a5a1317
# ╠═3779b846-e5ec-4239-a1d4-af2f8c2f10eb
# ╠═93df89f0-8429-4fcc-bd01-6982417f5134
# ╠═56afefe8-4452-4b2a-8a3b-e493ee1dd6c6
# ╟─fe9b889d-79c2-493b-9426-e33e6820cd90
# ╟─950d3b4e-f957-45b6-aa80-e3dfc765aad0
# ╟─50334069-a50c-467c-94ae-63b9b2264a18
# ╟─ddc141ba-d2e8-4ac4-8bc3-12fb1bb9fd4d
# ╟─4ad16c5c-73bc-4e42-9bfc-aea73a6bfbfe
# ╟─e89794b1-5bcd-4b6c-9cb2-77deca569c2e
# ╟─dcdb88ef-f04f-4ee8-87cc-bb26f396f064
# ╠═d9f5de8a-f3a2-41c9-9f3c-a0c8347368a4
# ╠═e6589d31-dce7-42c3-b494-db03fe561ae9
# ╟─7dbe9349-8b9e-4ac2-b4bf-b59f58a10ebc
# ╠═cf9e03db-42b7-41f6-80ce-4b12ddb93211
# ╠═1f27f695-9d83-4744-b6f0-e3280b746691
# ╠═bc90290e-81ba-4963-8a1f-5be15c0f3940
# ╠═9e847b90-54ae-4afb-bf72-1b9adeb97abc
# ╠═066b7505-e21b-467e-86c1-cea1ff80246e
# ╟─3f8db202-ac50-462d-b96d-ba629ca43325
# ╟─f7632708-7a0f-4c18-9f0b-44fb49aaeaa0
# ╟─cb45ae16-8893-4b80-bae0-fcef646a56b2
# ╟─99758072-3774-4bdc-9730-f0efff0e9330
# ╟─95cfdfe4-7fd8-4155-8a29-5a9711ab0f10
# ╟─f0eb5d23-e4c5-4189-a4aa-adf3219227cf
# ╟─b9b223de-0ff8-436e-9a4c-e056e1c3a412
# ╟─12866252-a5c6-43d0-92f1-d52df5a2d949
# ╟─a141cd27-6ea0-4f73-80b5-72d8e5770ed4
# ╟─d80f94c4-03d2-4aac-90f5-9415405b4412
# ╟─043f7a23-3b59-4e34-a8d3-9853cc66c228
# ╟─1596bc9f-f7e4-4d3d-9978-9da4eecbaede
# ╟─675d0bb0-4601-4f4e-bc7d-5d5fb2d70b18
# ╟─afea37f1-70c2-4aae-94f6-34cf7c1d9f8e
# ╟─ad8edd69-c164-4221-bdee-e7c9381ffcab
# ╟─9625798a-67df-49e4-91ce-c7e23ed2a177
# ╠═e4100b5b-b255-48db-a989-016fa72f8da5
# ╠═ef14041a-9753-456d-a242-d08dc0328507
# ╠═ecdc5f61-6041-42ef-819c-1d83c062c8e3
# ╠═7f6b0fe8-8ae1-4a42-87ea-2d5d3ae95181
# ╠═41fa1014-0d0a-4b30-9452-c5a1fc8e58b5
# ╠═553bce96-a29f-4a77-a193-8005914f4bfa
# ╠═e5d3b132-baca-4cdb-aeb0-09272610ed6f
# ╟─240a8494-158f-4db2-9d0d-c141b50dcd5d
# ╟─48921b9a-d3fd-4d03-9170-a414797d8dee
# ╟─693397a9-3c80-41e3-b618-fa9d5ecb42b7
# ╠═d5a956d0-b443-4ee2-9045-14146012a435
# ╠═2012c352-3bff-4c84-bc2e-d83b22d95349
# ╠═6b819b32-3601-48cc-9ff4-1df3e88e8034
# ╟─0858af84-9c5c-43b4-b5c9-fec1f1c4dba4
# ╠═ab8a18b1-4975-4f4a-994b-b08773489baf
# ╠═79064f4a-f4ad-4837-898a-d7a7e03f04da
# ╠═522d5f7b-7774-40f3-aa66-6dba2b48783a
# ╠═31a3d715-4a80-4696-96c3-8e5152797f9f
# ╠═7a4fbc79-9221-472f-acfa-d32249c8e828
# ╠═59ed7708-3c09-431e-8a09-1e94b0b90e2a
# ╠═cf516855-40fa-44d0-8d63-c145ef4577d0
# ╠═6734cf5a-a3dc-4de3-9f33-45e149a27ebc
# ╠═d9a98759-e596-4c81-98f5-2d4efdc964f6
# ╠═20c7902d-ebbf-49d1-a103-fb7524489524
# ╟─8c14d4a9-e2e7-4b18-bbba-a290872c2ca6
# ╠═9db5401c-8501-45b5-8d8e-070afcdfd867
# ╟─bf7ea880-5dff-419d-ba69-b4271b04bbe2
# ╟─5dc80a12-9265-4e5d-92b3-be8883f65ee7
# ╟─def4bd87-532b-4552-8229-620a27c04641
# ╟─1b9fe2e0-6a17-4f6b-bb23-f371e261e86e
# ╠═3eadc95c-b063-4ca6-ab5f-20acf5c56c2e
# ╠═4886158c-8dca-4d70-9903-7287624983ef
# ╠═33a0fc85-1069-4355-b2a6-6165f91d24e7
# ╠═ca18c7c8-6c58-4bd5-a92b-ed044c012603
# ╠═bced5a4c-9e43-4910-b233-738675a0f887
