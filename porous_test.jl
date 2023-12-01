include("TidalLoveNumbersPorous.jl")
using .TidalLoveNumbers

#enceladus core model:




n = 2
G = 6.67408e-11
e = 0.0047
ρₛ = [2422, 2422] 
r = [0, 0.1, 252.1-38-23] .* 1e3
μ = [1, 1+0im] .* 1e9
κ = [10e9, 10e9]
η = [1.0, 1.0] .* 1e25


bp = 2
tp = 2


# To recover solid-body behaviour
# set ρₗ=0, ϕ=0, α=0

ρₗ = [0, 1e3] 
# α  = [0, 0, 0.95, 0., 0]
κₗ = [0, 2.2e9]
ω = 2π / (33*60*60)
ηₗ = [0, 1.9e-3]
ϕ =  [0, 0.2]


ρ = (1 .- ϕ) .* ρₛ + ϕ .* ρₗ # bulk density

r = expand_layers(r)

g = get_g(r, ρ)



Ic = get_Ic(r[end,1], ρ[1], g[end,1], μ[1], "solid", 8, 4)

# Bprod = get_B_product(r, ρ, g, μ, κ, η, ω)[:,:,end,end]
Bprod1 = get_B_product(r, ρ, g, μ, κ, η, ω, ρₗ, κₗ, ηₗ, ϕ)[:, :, :, :]
Bprod2 = get_B_product(r, ρ, g, μ, κ, η, ω, ρₗ, κₗ, ηₗ, ϕ, bp, tp)[:, :, :, :]

Bprod4 = get_B_product(r, ρ, g, μ, κ, η, ω, ρₗ, κₗ, ηₗ, ϕ, 2, tp)[:, :, :, :]

# display(Bprod1)

# # P3 = zeros(3, 8)
# # P3[1,3] = 1
# # P3[2,4] = 1
# # P3[3,6] = 1

P4 = zeros(4, 8)
P4[1,3] = 1
P4[2,4] = 1
P4[3,6] = 1
P4[4,8] = 1

Pl = zeros(Int64, 8,4)
Pl[7,4] = 1

# # display(P4*Bprod4*Ic)
# # display(P4*Bprod2*Pl)
# # # display(Bprod4)

# # Ps = zeros(3, 4)
# # Ps[1,1] = 1
# # Ps[2,2] = 1
# # Ps[3,4] = 1

# # Bp = zeros(8)
# # Bp[7] = 1.0

# # y1_T1 = Bprod1*Ic

# # # display(y1_T1)

# # y_tp = Bprod4*Ic #+ Bprod2*Bp
# # # display(P4*Bprod4*Ic)
# # # display(Bprod2)
# # # display(T1)

# # # display(P3*Bprod1*Ic)
# # # display((Bprod3*Bprod2)*Bp)
# # # display(Bprod3)

# # # Bbp = get_B(r[bp], r[bp+1], g[bp], g[bp+1], ρ[bp], μ[bp], κ[bp], η[bp], ω, ρₗ[bp], κₗ[bp], ηₗ[bp], ϕ[bp])
# # # display(P3*Bprod1*Ic)
# # # y1 = Bprod*Ic
# # # # display(y1)

# # P1 = zeros(Float64, 3, 6)
# # P1[1,3] = 1
# # P1[2,4] = 1
# # P1[3,6] = 1

# # b = zeros(Float64, 3)
# # b[3] = (2n+1)/r[end,end] 

# # C = (P3* y1_T1) \ b

# # y_tp = Bprod4*Ic #+ Bprod2*Bp
# # # display(Bprod1*Ic*C)
# # # display(Bprod2)
# # # # display(T1)


# # # C = (P1 * y1) \ b

# # # Seems to be equivalent to using P1 and P2 
# # # y = C[1]*y1[:,1] + C[2]*y1[:,2] + C[3]*y1[:,3]
# # y = C[1]*y1_T1[:,1] + C[2]*y1_T1[:,2] + C[3]*y1_T1[:,3]
# # # display(y)
# # # P2 = zeros(Float64, 3, 6)
# # # P2[1,1] = 1
# # # P2[2,2] = 1
# # # P2[3,5] = 1

# # # y1 = P1 * Bprod * Ic * C 
# # # y2 = P2 * Bprod * Ic * C 

# # # println("k2 = ", y[5] - 1)  
# # # println("h2 = ", -g[end]*y[1] )



# ############################################ Marc's method ###############################################

# # println( -g[end]*y[2]  )

# # Bprod = get_B_product(r, ρ, g, μ, κ, η, ω)[:,:,end,end]

# # display(Bprod)

# yR = Bprod1[:,:,end,end]*Ic
# ytp = Bprod4[:,:,end,tp]*Ic + Bprod2[:,:,end,tp]*Pl
# yR = Bprod4[:,:,end,tp]*Ic + Bprod2[:,:,end,tp]*Pl
yR = Bprod1[:,:,end,end]*(Ic + Pl)


# # display(yR)
display(P4*yR)


# C1v = zeros(4)
# C2v = zeros(4)
# C3v = zeros(4)
# C4v = zeros(4)

# C1v[1] = 1
# C2v[2] = 1
# C3v[3] = 1
# C4v[4] = 1

# # display(yR)
# # display(yR *C1v)
# # display(yR *C2v)
# # display(yR *C3v)
# # display(yR *C4v)

# y1R = yR * C1v
# y2R = yR * C2v
# y3R = yR * C3v
# y4R = yR * C4v 

# y1tp = ytp * C1v
# y2tp = ytp * C2v
# y3tp = ytp * C3v
# y4tp = ytp * C4v 



# # Construct boundary condition matrix
# M = zeros(ComplexF64, 4,4)

# # Row 1 - Radial Stress
# M[1, 1] = y1R[3]
# M[1, 2] = y2R[3]
# M[1, 3] = y3R[3]
# M[1, 4] = y4R[3]

# # Row 2 - Tangential Stress
# M[2, 1] = y1R[4]
# M[2, 2] = y2R[4]
# M[2, 3] = y3R[4]
# M[2, 4] = y4R[4]

# # Row 3 - Potential Stress
# M[3, 1] = y1R[6]
# M[3, 2] = y2R[6]
# M[3, 3] = y3R[6]
# M[3, 4] = y4R[6]

# #  Row 4 - Darcy flux (r = r_tp)
# M[4, 1] = y1tp[8]
# M[4, 2] = y2tp[8]
# M[4, 3] = y3tp[8]
# M[4, 4] = y4tp[8]

M = P4*yR

display(M)
# # M[:,]


# # display(ytp*C1v)

# # display(P4*yR*C1v)
# # M[1,1] = (P4*yR*C1v)[1]
# # M[1:end-1, 1] .= (P4*yR*C1v)[1:end-1]
# # M[1:end-1, 2] .= (P4*yR*C2v)[1:end-1]
# # M[1:end-1, 3] .= (P4*yR*C3v)[1:end-1]
# # M[1:end-1, 4] .= (P4*yR*C4v)[1:end-1]
# # M[end, 1] = (ytp*C1v)[8]
# # M[end, 2] = (ytp*C2v)[8]
# # M[end, 3] = (ytp*C3v)[8]
# # M[end, 4] = (ytp*C4v)[8]

# # M[1:end-1,:] = P4*yR
# # M[end, :] = ytp

# # display(M)

# # M[:, 3] = P4*y1*C3v

# # display(M)

# # # M[1,1] = (y1*C1v)[3]
# # # M[1,2] = (y1*C2v)[3]

# # # display(M)P4*
# # P4*
# # P4*
# # P4*
b = zeros(ComplexF64, 4)
b[3] = (2n+1)/r[end,end] 
C2 = M \ b

# # C3 = M[1:3,1:3] \ b[1:3]

y = yR*C2

# # display(ytp[:,1:3]*C3)
# # b[1] = (ytp[:,1:3]*C3)[3]
# # b[2] = (ytp[:,1:3]*C3)[4]
# # b[3] = (ytp[:,1:3]*C3)[6]


# # C4 = M \ b
display(y)
display(C2)
# # # display(C3)
# # # display(C4)

# # y1 = ytp*C2
# # y2 = ytp*C4

# # # y = C2[1]*yR[:,1] + C2[2]*yR[:,2] + C2[3]*yR[:,3]+ C2[4]*yR[:,4]

# # # y = C2[1]*ytp[:,1] + C2[2]*ytp[:,2] + C2[3]*ytp[:,3]+ C2[4]*ytp[:,4]
# # # display(C2)
# # display(y1)
# # display(y2)

# # # display(C)
# # display(C2)

k2 = y[5] - 1
println("k2 = ", k2)  
println("h2 = ", -g[end]*y[1] )

Ediss = -21/2 * imag(k2) * (ω*r[end,end])^5/G * e^2

println(Ediss/1e9 )





