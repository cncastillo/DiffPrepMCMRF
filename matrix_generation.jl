## Aux functions
""""Calculates the normalized moments Mₖ = 1/tᵏ ∫ᵗG(τ)τᵏ dτ at the end of the sequence. """
function get_Mmatrix(seq::Sequence; axis=1, τ_sample=dur(seq))
    τ = τ_sample # Seq Duration [ms]
    T0 = get_block_start_times(seq)
    M0, M1, M2, M3 = Float64[], Float64[], Float64[], Float64[]
    for i = 1:length(seq)
        #Gradient
        Gi = seq[i].GR[axis]
        N = length(Gi.A)
        delay = Gi.delay #Durations of delay [s]
        #Timings
        if N > 1
            δ = ones(N) * Gi.T / (N-1) #Durations of pulse [s]
            T = [sum(δ[1:j]) for j = 1:N-1]
            T = T0[i] .+ delay .+ [0; T] #Position of pulse
            #Moment calculations - P1 model
            append!(M0, δ/τ^1)
            append!(M1, δ.*(T)/τ^2)
            append!(M2, δ.*(T.^2 .+ δ.^2/6)/τ^3)
            append!(M3, δ.*(T.^3 .+ T .* δ.^2/2)/τ^4)
        end
    end
    [M0'; M1'; M2'; M3']
end

"""Slew rate matrix: |SR*g| ≤ Smax."""
function get_SRmatrix(seq::Sequence; axis = 1)
    SR = Bidiagonal{Float64}[]
    for i = 1:length(seq)
        #Gradient
        Gi = seq[i].GR[axis]
        N = length(Gi.A)
        if N > 1
            Δt = ones(N) * Gi.T / (N-1)
            dv = Δt
            ev = Δt[1:end-1]
            SRi = Bidiagonal(-1 ./ dv, 1 ./ ev, :U)
            push!(SR, SRi)
        end
    end
    SR
end

"""Maxwell matrix MX=∫g²dt."""
function get_MXmatrix(seq::Sequence; axis = 1)
    MX = Tridiagonal{Float64}[]
    τ = dur(seq) # Seq Duration [ms]
    for i = 1:length(seq)
        #Gradient
        Gi = seq[i].GR[axis]
        N = length(Gi.A)
        if N > 1
            Δt = ones(N) * Gi.T / (N-1)
            dd = Δt[2:end]/6
            d = 2Δt/3
            MXi = Tridiagonal(dd, d, dd)
            push!(MX, MXi/τ^3)
        end
    end
    MX
end

"Calculates the `b`-matrix, such as `b`-value = g' B g [s/mm2] with g [T/m]."
get_Bmatrix(seq::Sequence; axis=1) = begin
    T0 = get_block_start_times(seq)[1:end-1]
    #Calculating timings
    T = Float64[]
    δ = Float64[]
    for (i, s) = enumerate(seq)
        #Gradient
        Gi = s.GR[axis]
        N = length(Gi.A)
        delay = Gi.delay
        if N > 1
            δi = ones(N) * Gi.T / (N-1)
            Ti = [sum(δi[1:j]) for j = 1:N-1]
            Ti = T0[i] .+ delay .+ [0; Ti] #Position of pulse
            append!(δ, δi)
            append!(T, Ti)
        end
    end
    τ = dur(seq) + δ[end]
    Nsamples = length(T)
	ij = [max(i,j) for i=1:Nsamples, j=1:Nsamples]
	α = [(i==j) ? 2/3 : 1/2 for i=1:Nsamples, j=1:Nsamples]
	b = (δ' .* δ) .* (τ .- T[ij] .- α .* δ[ij])
	b_value = (2π*γ)^2*1e-6*b # Trace of B tensor
	b_value
end

## TO SCANNER (Philips)
"""Duration in [s] => samples, with dwell-time of Δt = 6.4 μs."""
δ2N(δ) = floor(Int64, δ * 156250) + 2

"""Exports diffusion preparation waveforms for their use in the scanner."""
function write_diffprep_fwf(G1, G2, G3, bmax, Gmax, Smax; filename="./qte_vectors_input.txt", name="Maxwell2",
    precision::Int=6, dwell_time=6.4e-6, verbose=false)
    open(filename, "w") do io
        t1 = range(0, G1.GR.dur[1] - maximum(G1.GR.delay), step=dwell_time) #length=δ2N(maximum(G1.GR.T))) #step=dwell_time) #
		t2 = range(0, G2.GR.dur[1] - maximum(G2.GR.delay), step=dwell_time) #length=δ2N(maximum(G2.GR.T)))
        t3 = range(0, G3.GR.dur[1] - maximum(G3.GR.delay), step=dwell_time) #length=δ2N(maximum(G3.GR.T)))
        maxN = max(length(t1), length(t2), length(t3))
        Gx1, Gy1, Gz1 = KomaMRIBase.get_grads(G1, Array(t1).+maximum(G1.GR.delay))
		Gx2, Gy2, Gz2 = KomaMRIBase.get_grads(G2, Array(t2).+maximum(G2.GR.delay))
        Gx3, Gy3, Gz3 = KomaMRIBase.get_grads(G3, Array(t3).+maximum(G3.GR.delay))
        Gx1_round = round.(Gx1 ./ Gmax, digits=precision)
        Gx2_round = round.(Gx2 ./ Gmax, digits=precision)
        Gx3_round = round.(Gx3 ./ Gmax, digits=precision)
        Gy1_round = round.(Gy1 ./ Gmax, digits=precision)
        Gy2_round = round.(Gy2 ./ Gmax, digits=precision)
        Gy3_round = round.(Gy3 ./ Gmax, digits=precision)
        Gz1_round = round.(Gz1 ./ Gmax, digits=precision)
        Gz2_round = round.(Gz2 ./ Gmax, digits=precision)
        Gz3_round = round.(Gz3 ./ Gmax, digits=precision)
        if verbose
        println("Δt1=$(t1[2]-t1[1]) $(Gx1_round[1]) $(Gx1_round[end]) $(Gy1_round[1]) $(Gy1_round[end]) $(Gz1_round[1]) $(Gz1_round[end])")
        println("Δt2=$(t2[2]-t2[1]) $(Gx2_round[1]) $(Gx2_round[end]) $(Gy2_round[1]) $(Gy2_round[end]) $(Gz2_round[1]) $(Gz2_round[end])")
        println("Δt3=$(t3[2]-t3[1]) $(Gx3_round[1]) $(Gx3_round[end]) $(Gy3_round[1]) $(Gy3_round[end]) $(Gz3_round[1]) $(Gz3_round[end])")
        end
        M01 =  [sum(floor.(Int32, Gx1_round*10^precision)) sum(floor.(Int32, Gy1_round*10^precision)) sum(floor.(Int32, Gz1_round*10^precision))]
        M02 = -[sum(floor.(Int32, Gx2_round*10^precision)) sum(floor.(Int32, Gy2_round*10^precision)) sum(floor.(Int32, Gz2_round*10^precision))]
        M03 =  [sum(floor.(Int32, Gx3_round*10^precision)) sum(floor.(Int32, Gy3_round*10^precision)) sum(floor.(Int32, Gz3_round*10^precision))]
        M0 = M01 .+ M02 .+ M03
        Gx1_diff =  maximum(abs.(Gx1_round[2:end] .- Gx1_round[1:end-1]))
        Gx2_diff =  maximum(abs.(Gx2_round[2:end] .- Gx2_round[1:end-1]))
        Gx3_diff =  maximum(abs.(Gx3_round[2:end] .- Gx3_round[1:end-1]))
        Gy1_diff =  maximum(abs.(Gy1_round[2:end] .- Gy1_round[1:end-1]))
        Gy2_diff =  maximum(abs.(Gy2_round[2:end] .- Gy2_round[1:end-1]))
        Gy3_diff =  maximum(abs.(Gy3_round[2:end] .- Gy3_round[1:end-1]))
        Gz1_diff =  maximum(abs.(Gz1_round[2:end] .- Gz1_round[1:end-1]))
        Gz2_diff =  maximum(abs.(Gz2_round[2:end] .- Gz2_round[1:end-1]))
        Gz3_diff =  maximum(abs.(Gz3_round[2:end] .- Gz3_round[1:end-1]))
        SRx1 = Gx1_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRx2 = Gx2_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRx3 = Gx3_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRy1 = Gy1_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRy2 = Gy2_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRy3 = Gy3_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRz1 = Gz1_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRz2 = Gz2_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        SRz3 = Gz3_diff * Gmax / dwell_time # 6.4 us is the dwell-time
        if verbose
        println("SR1 = [$SRx1, $SRy1, $SRz1]")
        println("SR2 = [$SRx2, $SRy2, $SRz2]")
        println("SR3 = [$SRx3, $SRy3, $SRz3]")
        end
        println("SR = [$(max(SRx1, SRy1, SRz1)), $(max(SRx2, SRy2, SRz2)), $(max(SRx3, SRy3, SRz3))]")
        @assert (SRx1 <= Smax) && (SRx2 <= Smax) && (SRx3 <= Smax)
        @assert (SRy1 <= Smax) && (SRy2 <= Smax) && (SRy3 <= Smax)
        @assert (SRz1 <= Smax) && (SRz2 <= Smax) && (SRz3 <= Smax)
        if verbose
        println("M01 = [$(M01[1]), $(M01[2]), $(M01[3])]")
        println("M02 = [$(M02[1]), $(M02[2]), $(M02[3])]")
        println("M03 = [$(M03[1]), $(M03[2]), $(M03[3])]")
        end
        println("M0 = $(M0.*10.0^(-precision))")
        # @assert all(M0 .== 0)
        MX1 =  [sum(floor.(Int32, Gx1_round*10^precision).^2) sum(floor.(Int32, Gy1_round*10^precision).^2) sum(floor.(Int32, Gz1_round*10^precision).^2)]
        MX2 = -[sum(floor.(Int32, Gx2_round*10^precision).^2) sum(floor.(Int32, Gy2_round*10^precision).^2) sum(floor.(Int32, Gz2_round*10^precision).^2)]
        MX3 =  [sum(floor.(Int32, Gx3_round*10^precision).^2) sum(floor.(Int32, Gy3_round*10^precision).^2) sum(floor.(Int32, Gz3_round*10^precision).^2)]
        MX = MX1 .+ MX2 .+ MX3
        if verbose
        println("MX1 = [$(MX1[1]), $(MX1[2]), $(MX1[3])]")
        println("MX2 = [$(MX2[1]), $(MX2[2]), $(MX2[3])]")
        println("MX3 = [$(MX3[1]), $(MX3[2]), $(MX3[3])]")
        end
        println("MX = $(abs.(MX).*10.0^(-2precision))")
        # @assert all(MX .== 0)
        #BVAL
        t = range(0, dur(G1+G2+G3), step=dwell_time)
        Gx, Gy, Gz = KomaMRIBase.get_grads(G1-G2+G3, Array(t))
        bvalx = (2π*γ)^2 * 1e-6 * sum(cumsum(Gx * dwell_time).^2 * dwell_time)
        bvaly = (2π*γ)^2 * 1e-6 * sum(cumsum(Gy * dwell_time).^2 * dwell_time)
        bvalz = (2π*γ)^2 * 1e-6 * sum(cumsum(Gz * dwell_time).^2 * dwell_time)
        bval = round(bvalx+bvaly+bvalz, digits=3)
        println("bval_calc = [$bvalx $bvaly $bvalz] ($bval s/mm2)")
        #Header
        N1, N2, N3 = length(t1), length(t2), length(t3)
        println("N1 = $N1 N2 = $N2 N3 = $N3")
        date = "#Generated on $(now())\r\n"
        vars =  @sprintf "%s %s %s %s %s %s %s\r\n" "#Name"*" "^(length(name)-5) "N1"*" "^(length(string(N1))-1) "N2"*" "^(length(string(N2))-1) "N3"*" "^(length(string(N3))-1) "bval"*" "^(length(string(round(bmax,digits=1)))-3) "Gmax"*" "^(length(string(round(Gmax,digits=1)))-3) "Smax"
        unit =  @sprintf "%s %s %s %s\r\n" "#"*" "^(length(name)+length(string(N1))+length(string(N2))+length(string(N3))+2)  "s/mm2"*" "^(length(string(bval))-3) "mT/m"*" "^(length(string(round(Gmax,digits=1)))-3) "T/m/s"
        line =  @sprintf "%s %i %i %i %.1f %.1f %.1f\r\n" name N1 N2 N3 bval Gmax*1e3 Smax
        write(io, date)
        write(io, vars)
        write(io, unit)
        write(io, line)
        for i = 1:maxN
            fx1, fy1, fz1 = i ≤ length(t1) ? (Gx1_round[i], Gy1_round[i], Gz1_round[i]) : (0,0,0)
            fx2, fy2, fz2 = i ≤ length(t2) ? (Gx2_round[i], Gy2_round[i], Gz2_round[i]) : (0,0,0)
            fx3, fy3, fz3 = i ≤ length(t3) ? (Gx3_round[i], Gy3_round[i], Gz3_round[i]) : (0,0,0)
            line = @sprintf "% .6f % .6f % .6f % .6f % .6f % .6f % .6f % .6f % .6f\r\n" fx1 fy1 fz1 fx2 fy2 fz2 fx3 fy3 fz3
            write(io, line)
        end
    end
end