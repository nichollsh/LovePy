include("TidalLoveNumbers.jl")
using .TidalLoveNumbers

#enceladus test model:
n = 2
# ρ = [5376.81502017327, 3300, 1000, 1000, 940]
# r = [0, 720, 1441, 1541, 1556, 1561] .* 1e3
# μ = [40, 40, 1.0, 3.5, 3.5+0im] .* 1e9
# κ = [5000e9, 5000e9, 5000e9, 5000e9, 5000e9] * 1e8
# η = [1.0, 1.0, 1.0, 1.0, 1.0] .* 1e24


# ρₗ = [0, 0, 3000, 0, 0]
# α  = [0, 0, 0.1, 0., 0]
# κₗ = [0, 0, 1e10, 0, 0]
# ω = 2.047e-5
# ηₗ = [0, 0, 1e-2, 0, 0]
# ϕ = [0, 0, 0.0, 0, 0]


G = 6.67408e-11
e = 0.0047
ρₛ = [2422, 2422.0]
r = [0, 0.01, 252.1-38-23] .* 1e3
μ = [1, 1+0im] .* 1e9
κ = [10e9, 10e9]
η = [1.0, 1.0] .* 1e25
ω = 2π / (33*60*60)



# r, ρ, μ, κ, η = expand_layers(r, ρ, μ, κ, η)




r = expand_layers(r)





# View properties of high res layers
# display(hcat(r, ρ, μ, κ))



g = get_g(r, ρₛ)

# display(r)
# display(g)

Ic = get_Ic(r[end,1], ρₛ[1], g[end,1], μ[1], "solid")


Bprod = get_B_product(r, ρₛ, g, μ, κ, η, ω)[:,:,end,end]

# display(Bprod)

y1 = Bprod*Ic



P1 = zeros(Float64, 3, 6)
P1[1,3] = 1
P1[2,4] = 1
P1[3,6] = 1

b = zeros(Float64, 3)
b[3] = (2n+1)/r[end] 

display(P1*y1)

C = (P1 * y1) \ b

# Seems to be equivalent to using P1 and P2 
y = C[1]*y1[:,1] + C[2]*y1[:,2] + C[3]*y1[:,3]

# display(P1*y)

P2 = zeros(Float64, 3, 6)
P2[1,1] = 1
P2[2,2] = 1
P2[3,5] = 1


# y1 = P1 * Bprod * Ic * C 
# y2 = P2 * Bprod * Ic * C 


# display(y1)
# display(y2)
# display(b)


println("k2 = ", y[5] - 1)  
println("h2 = ", -g[end]*y[1] )

############################################ Marc's method ###############################################

# println( -g[end]*y[2]  )

# Bprod = get_B_product(r, ρ, g, μ, κ, η, ω)[:,:,end,end]

# # display(Bprod)

# y1 = Bprod*Ic

C1v = zeros(3)
C2v = zeros(3)
C3v = zeros(3)

C1v[1] = 1
C2v[2] = 1
C3v[3] = 1

M = zeros(ComplexF64, 3,3)
M = P1*y1
# M[:, 1] = P1*y1*C1v
# M[:, 2] = P1*y1*C2v
# M[:, 3] = P1*y1*C3v

# display(M)

# display(P1*y1)

# # M = zeros(ComplexF64, 3,3)
# # M[1, :] = P1*y1*C1v
# # M[:, 2] = P1*y1*C2v
# # M[:, 3] = P1*y1*C3v


# # M = zeros(Comp
# # M[1,2] = (y1*C2v)[3]

# # display(M)

# C2 = M \ b

# # display(y)
# # y = C2[1]*y1[:,1] + C2[2]*y1[:,2] + C2[3]*y1[:,3]
# display(y)
# # display(C)
# display(C2)





