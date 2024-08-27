# Code used to generate moment-compensated diffusion gradient waveforms
# Sequence optimization for diffusion prepared motion-compensated MRF

# Ease of use, and plots
using KomaMRI 
# Optimization
using JuMP, Ipopt
# Matrices
using LinearAlgebra: I, Bidiagonal, norm, Diagonal, Tridiagonal
# To save results    
using Printf, Dates

include("matrix_generation.jl")

########################################################################
# Optimization setup
########################################################################

# Parameters
dwell_time = 6.4e-6 # s
Gmax = 62e-3        # mT/m
Smax = 70           # mT/m/ms
k = 1               # moment nulling
n_dwells = 4        # Node separation
maxwell = true      # Maxwell/concomitant gradient compensation
gap_left_ms = 1     # ms, gap to the left of the RF
gap_right_ms = 0    # ms, gap to the right of the RF
# RF pulse timings
adia = "HS2"        # Adiabatic pulse name
δ1 = 7.2064e-3      # s, First gradient duration
δ2 = 14.9952e-3     # s, Second gradient duration
δ3 = 7.2064e-3      # s, Third gradient duration
Δ1 = 17.5024e-3     # s, Time between start of grad1 and start of grad2
Δ2 = 42.5076e-3     # s, Time between start of grad1 and start of grad3

########################################################################
# Ensure correct sampling of gradient waveforms
########################################################################

# Match gradient timings with dwell time
gap_left  = floor(Int64, gap_left_ms * 1e-3 / (n_dwells * dwell_time))
gap_right = floor(Int64, gap_right_ms * 1e-3 / (n_dwells * dwell_time))
δ1_new = floor(Int64, δ1 / dwell_time) * dwell_time # Making the waveform match the dwell time
δ2_new = floor(Int64, δ2 / dwell_time) * dwell_time # Making the waveform match the dwell time
δ3_new = floor(Int64, δ3 / dwell_time) * dwell_time # Making the waveform match the dwell time
@assert δ1_new ≈ δ1 "δ1_new = $(δ1_new*1e3) != δ1 = $(δ1*1e3)"
@assert δ2_new ≈ δ2 "δ2_new = $(δ2_new*1e3) != δ2 = $(δ2*1e3)"
@assert δ3_new ≈ δ3 "δ3_new = $(δ3_new*1e3) != δ3 = $(δ3*1e3)"
rf1 = Δ1 - δ1
rf2 = Δ2 - δ2 - Δ1
# Calculating number of samples for each gradient
N1 = floor(Int64, δ1 / (n_dwells * dwell_time)) + 1; println("Nnodes_grad1 = $N1")
N2 = floor(Int64, δ2 / (n_dwells * dwell_time)) + 1 # δ1/N1 = δ2/N2
N3 = floor(Int64, δ3 / (n_dwells * dwell_time)) + 1
# Using an slew rate slightly below the maximum
Smax_discrete = Smax * 0.999
# Calculate max slew rate due to node distances 
dt = max(δ1 / (N1-1), δ2 / (N2-1), δ3 / (N3-1))
Smax_discrete = Gmax / (dt * ceil(Int64, (Gmax / Smax) / dt))
println("Smax_discrete = ", Smax_discrete)

########################################################################
# Generating optimization matrices 
########################################################################

# Generating sequence
DIF = Sequence()
global DIF =  Sequence([Grad(x -> 1e-3, δ1, N1; delay=0)])
global DIF += Sequence([Grad(x -> 1e-3, δ2, N2; delay=rf1)])
global DIF += Sequence([Grad(x -> 1e-3, δ3, N3; delay=rf2)])

# Optimization matrices (piecewise-linear to nodes)
B =  get_Bmatrix(DIF)           # B-value
SR = get_SRmatrix(DIF)          # Slew-rate matrices
M =  get_Mmatrix(DIF)[1:k+1,:]  # k Moment nulling
MX = get_MXmatrix(DIF)          # Maxwell matrices

########################################################################
# Solve optimization 
########################################################################

# Sequence name
seq_name = ""
seq_name = maxwell ? "$(seq_name)MX_MC$(k)" : "$(seq_name)MC$(k)"  #Name of the sequnce
seq_name = adia != "" ? "$(adia)_$(seq_name)" : seq_name       #Name of the sequnce
seq_name *= gap_left_ms > 0 ? "_gap" : ""
println("#################### $seq_name ####################")
# Optimization problem
model = Model(Ipopt.Optimizer)
set_silent(model)
# Variables
@variable(model, -Gmax <= g1[1:N1] <= Gmax, start=Gmax);
@variable(model, -Gmax <= g2[1:N2] <= Gmax, start=Gmax);
@variable(model, -Gmax <= g3[1:N3] <= Gmax, start=Gmax);
# Objective function
@objective(model, Max, [g1;-g2;g3]'*B*[g1;-g2;g3]);
# Constraints
@constraint(model, moments_final, M *[g1;-g2;g3] .== 0);
@constraint(model, slewrate, -Smax_discrete .<= [SR[1]*g1; -SR[2]*g2; SR[3]*g3] .<= Smax_discrete);
@constraint(model, ends, [g1[1]; g2[1:1+gap_right]; g3[1]; g1[N1-gap_left:N1]; g2[N2-gap_left:N2]; g3[N3]] .== 0)
if maxwell
    @constraint(model, concomitant, g1'*MX[1]*g1 - g2'*MX[2]*g2 + g3'*MX[3]*g3 == 0);
end
# Solve
optimize!(model)
# Retrieving solution
gx1 = value.(g1)
gx2 = value.(g2)
gx3 = value.(g3)
global gx = [gx1; -gx2; gx3]
# Optimization status
if termination_status(model) == MOI.LOCALLY_SOLVED
    println( "Solved! 😃" )
else
    println( "NOT Solved 😢" )
end
# Putting solution inside Sequence object for plotting
# Creating the RFs
t = range(-1.5*rf1/2, 1.5*rf1/2, 80)
β = 4e2 #frequency modulation param (rad/s)
B1 = 2*13.5e-6 * sech.(β * t)
R1 = [RF(B1, rf1, 0, δ1);;]
R2 = [RF(B1, rf2, 0, δ2);;]
# Sequence with solution
global DIF =  Sequence([Grad(gx1,δ1,0); Grad(gx1,δ1,0); Grad(gx1,δ1,0);;],R1)
global DIF += Sequence([Grad(gx2,δ2,0); Grad(gx2,δ2,0); Grad(gx2,δ2,0);;],R2)
global DIF += Sequence([Grad(gx3,δ3,0); Grad(gx3,δ3,0); Grad(gx3,δ3,0);;])
bmax = 3 * objective_value(model)

########################################################################
# Display results
########################################################################

# Flip the gradient if the first gradient's x area is negative
inv = sum(DIF[1].GR[1].A) <= 0
global DIF = inv ? -DIF : DIF
# Plots
τ = dur(DIF) * 1e3
p1 = plot_seq(DIF; slider=false, range=[0,τ], title="$seq_name $(round(bmax, digits=2)) s/mm2")
p2 = plot_M0(DIF;  slider=false, range=[0,τ])
p3 = plot_M1(DIF;  slider=false, range=[0,τ])
p4 = plot_slew_rate(DIF; slider=false, range=[0,τ])
p = [p1; p2; p3; p4]
display(p)
# Save resulting waveform to TXT
write_diffprep_fwf(DIF[1], DIF[2], DIF[3], bmax, Gmax, Smax;
        filename="$seq_name.txt", name=seq_name, verbose=false)
println( "λ0 = $(abs(round(M[1,:]'*gx/Gmax,digits=3))), λ1 = $(abs(round(M[2,:]'*gx/Gmax,digits=3)))" )
println( "MX ∫g1²-∫g2²+∫g3²=$(gx1'*MX[1]*gx1 - gx2'*MX[2]*gx2 + gx3'*MX[3]*gx3)")
println( "b-value: $(round(bmax, digits=2)) s/mm2" )
println( seq_name )
println("Finished! 💃")