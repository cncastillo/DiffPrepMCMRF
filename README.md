# Optimization of moment-compensated diffusion gradient waveforms

![](DiffPrep_zoom.svg)

## Install
Download the repo:
```bash
git clone git@github.com:cncastillo/DiffPrepMCMRF.git
cd DiffPrepMCMRF
```
## Run
To run the waveform optimization, use the following command inside of the Julia REPL (`julia`):
```julia
julia> include("DiffPrepMCMRF.jl")
```

## Expected output
```julia
Nnodes_grad1 = 282
Smax_discrete = 69.07352194874397
#################### HS2_MX_MC1_gap ####################

******************************************************************************
This program contains Ipopt, a library for large-scale nonlinear optimization.
 Ipopt is released as open source code under the Eclipse Public License (EPL).
         For more information visit https://github.com/coin-or/Ipopt
******************************************************************************

Solved! 😃
[ Info: Listening on: 127.0.0.1:8492, thread id: 1
SR = [69.08156250000012, 69.08156250000012, 69.08156250000107]
M0 = [-0.001671 -0.001671 -0.001671]
MX = [0.069523650485 0.069523650485 0.069523650485]
bval_calc = [283.72156924233957 283.72156924233957 283.72156924233957] (851.165 s/mm2)
N1 = 1127 N2 = 2344 N3 = 1127
λ0 = 0.0, λ1 = 0.0
MX ∫g1²-∫g2²+∫g3²=-1.3074541449498156e-12
b-value: 851.17 s/mm2
HS2_MX_MC1_gap
Finished! 💃
```
Diffusion prepared module:
![image](https://github.com/user-attachments/assets/ef06b817-ecef-4572-b33d-53ad92a76bb2)
Waveforms stored as TXT:
![image](https://github.com/user-attachments/assets/76b62ca9-b05b-4be0-9e5e-1082ffd3263b)


