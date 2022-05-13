using Gridap
using GridapGmsh
using Gridap.Fields
using Gridap.CellData
using Gridap.Arrays
using LineSearches: BackTracking, Static, MoreThuente
using FillArrays
using Gridap.Fields: meas
using Statistics
using JLD2

"""
Jet analysis
Inlet velocity
"""



Re = 1000
D = 0.1 #[m] diameter
ν = 0.001 # [m2/s] kinematic viscosity


u0 = Re*ν/D

order = 1 #Order of pressure and velocity

hf = VectorValue(0.0,0.0) #external force


#ODE settings
t0 = 0.0
dt = 0.01
tF = 3

Ntimestep = (tF-t0)/dt
initial_condition = false #print model of initial condition



#MESH DEFINITION

model = GmshDiscreteModel("Jet1.msh")

writevtk(model,"model")

#Function of velocity and pressure on the boundaries
u_wall(x, t) = VectorValue(0,0)
u_wall(t::Real) = x -> u_wall(x, t)

u_inlet(x, t) = VectorValue(u0, 0) #inlet velocity Function
u_inlet(t::Real) = x -> u_inlet(x, t)

p0(x, t) = 0
p0(t::Real) = x -> p0(x, t)


u_start(x, t) = (abs(x[2])>D/2) ? VectorValue(0, 0) : VectorValue(u0, 0)
u_start(t::Real) = x -> u_start(x, t)


reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=["P_inlet","Inlet","Walls","P_walls"]) #diri0 union of the walls, diri1 top wall
reffeₚ = ReferenceFE(lagrangian,Float64, order)
Q = TestFESpace(model,reffeₚ, conformity=:H1, dirichlet_tags=["P_outlet","Outlet"])

U = TransientTrialFESpace(V, [u_inlet, u_inlet, u_wall, u_wall])
P = TransientTrialFESpace(Q, [p0,p0]) 



Y = MultiFieldFESpace([V, Q])
X = TransientMultiFieldFESpace([U, P])

degree = 4*order
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)


h = lazy_map(h->h^(1/2),get_cell_measure(Ω)) #get the area of each cell (because in 2D), then ^1/2



# Momentum residual, without the viscous term
Rm(t,(u,p)) = ∂t(u) + u⋅∇(u) + ∇(p) - hf

# Continuity residual
Rc(u) = ∇⋅u


function τ(u,h)
   
    β=u0 
    τ₂ = h^2/(4*ν)
    val(x) = x
    val(x::Gridap.Fields.ForwardDiff.Dual) = x.value
    u = val(norm(u))
    
    if iszero(u)
        return τ₂
        
    end
    τ₃ =  dt/2 #h/(2*u) #0  dt/2 #

    τ₁ = h/(2*u) #h/(2*u) #
    return 1/(1/τ₁ + 1/τ₂ + 1/τ₃)
    
end


#τb(u,h) = (u⋅u)*τ(u,h)
τb(u,h) = (u⋅u)*τ(u,h)

var_equations(t,(u,p),(v,q)) = ∫(
    ν*∇(v)⊙∇(u) # Viscous term
    + v⊙Rm(t,(u,p)) # Other momentum terms
    + q*Rc(u)
 )dΩ # Continuity


stab_equations(t,(u,p),(v,q)) = ∫(  (τ∘(u,h)*(u⋅∇(v) + ∇(q)))⊙Rm(t,(u,p)) # First term: SUPG, second term: PSPG
    +τb∘(u,h)*(∇⋅v)⊙Rc(u) # Bulk viscosity. Try commenting out both stabilization terms to see what happens in periodic and non-periodic cases
)dΩ


res(t,(u,p),(v,q)) = var_equations(t,(u,p),(v,q)) + stab_equations(t,(u,p),(v,q))
op = TransientFEOperator(res,X,Y)


nls = NLSolver(show_trace=true, method=:newton, linesearch=MoreThuente(), iterations=30)

solver = FESolver(nls)




U0 = U(0.0)
P0 = P(0.0)
X0 = X(0.0)

#initial condition
uh0 = interpolate_everywhere(u_start(0), U0)
ph0 = interpolate_everywhere(0, P0)
xh0 = interpolate_everywhere([uh0, ph0], X0)

#initial condition - derivative; It is not ideal, the first iteration are not perfect
vuh0 = interpolate_everywhere(VectorValue(0,0), U0)
vph0 = interpolate_everywhere(0, P0)
vxh0 = interpolate_everywhere([vuh0, vph0], X0)

writevtk(Ω, "initial_condition", cellfields=["uh" => uh0, "ph" => ph0])

#For extracting results, for fft
N_samples = 100
#U_vector = zeros(length(uh0.free_values), N_samples)

ρ∞ = 0.8 #ρ∞=1 no dissipation, ρ∞=0 max dissipation, ρ∞=0.5 quite good 
ode_solver = GeneralizedAlpha(nls,dt,ρ∞)
sol_t = solve(ode_solver,op,(xh0,vxh0),t0,tF)


_t_nn = t0
iteration = 0
s_iteration = 1

createpvd("Jet_2d") do pvd
  for (xh_tn, tn) in sol_t
    global _t_nn
    _t_nn += dt
    global iteration
    global s_iteration
    global U_vector
    iteration += 1
    println("it_num = $iteration\n")
    uh_tn = xh_tn[1]
    ph_tn = xh_tn[2]
    ωh_tn = ∇ × uh_tn
    """
    if tn>100
      U_vector[: , s_iteration] =  uh_tn.free_values
      s_iteration += 1
    end
    """
    if mod(iteration,1  )<1 #useful for printing not every iteration
      pvd[tn] = createvtk(Ω, "Results/Jet_2d_$_t_nn" * ".vtu", cellfields=["uh" => uh_tn, "ph" => ph_tn, "wh" => ωh_tn])
    end
    
  end

end

#@save "Jet.jld2" U_vector
