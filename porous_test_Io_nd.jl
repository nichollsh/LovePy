include("TidalLoveNumbersPorous.jl")
using .TidalLoveNumbers
using DoubleFloats
using PyPlot
using Statistics

prec = TidalLoveNumbers.prec
precc = TidalLoveNumbers.precc

non_dim = true

G = prec(6.6743e-11)
e = 0.0041

h_core = 700.0
h_mantle_low = 800. + 0.
h_mantle_up = 300. - 0.

h_crust = 20.
ω0 = 2*2.047e-5
ω = ω0

#enceladus test model:
n = 2
ρₛ = [7640, 3300, 3300, prec(3300)]
r = [0, 
     h_core, 
     h_core+h_mantle_low, 
     h_core+h_mantle_low+h_mantle_up, 
     h_core+h_mantle_low+h_mantle_up+h_crust] .* 1e3
μ = [40, 60, 60, prec(60)] .* 1e9
κ = [100e9, 100e9, 100e9, 100e9]
η = [1e23, 1e23, 1e23, 1e25]

R = r[end]
μ0 = μ[end]
bp = 3
tp = 3


ρₗ = [0, 0, 3300, 0]
# α  = [0, 0, 0.95, 0., 0]
κₗ = [0, 0, 10e10, 0]
k = [0, 0, 1e-6, 0]

ηₗ = [0, 0, 1e-2, 0]
ϕ =  [0, 0, prec(0.3), 0]

ρ = (1 .- ϕ) .* ρₛ + ϕ .* ρₗ # bulk density
r = expand_layers(r)
g = get_g(r, ρ)

 
ρ0 = ρ[end]
μ0 = ρ0*g[end,end]*R

Gnd = G*ρ0^2*R^2/μ0;



# Now non-dimensionalise?
T = 2π/ω0

if non_dim
    r = r./R
    ρₛ = ρₛ ./ ρ0
    ρ = ρ ./ ρ0
    μ = μ ./ μ0
    κ = κ ./ μ0
    η = η ./ (μ0 * 2π/ω0 )

    ρₗ = ρₗ ./ ρ0
    κₗ = κₗ ./ μ0
    k = k ./ R^2

    ηₗ = ηₗ ./ (μ0 * 2π/ω0 )


    ω = 2π

    set_G(Gnd)
    g = get_g(r, ρ)

end







ηs = 10 .^ collect(19:1:19)

Edot1 = zeros(length(ηs))
Edot2 = zeros(length(ηs))
k2_1 = zeros(ComplexF64, length(ηs))
k2_2 = zeros(ComplexF64, length(ηs))

yall = zeros(precc, 8, size(r)[1]-1,size(r)[2])
for i in eachindex(ηs)
    η[3] = ηs[i] 
    if non_dim
        η[3] = η[3]/(μ0 * T )
    end

    μc =  1im*ω*μ ./ (1im*ω .+ μ ./ η)

    Ic = get_Ic(r[end,1], ρ[1], g[end,1], μ[1], "liquid")

    Bprod1 = get_B_product(r, ρ, g, μc, κ)[:, :, :, :]

    # Projection matrix for the third, fourth, sixth, and eigth 
    # components of the solution
    P1 = zeros(3, 6)
    P1[1,3] = 1
    P1[2,4] = 1
    P1[3,6] = 1

    yR = Bprod1[:,:,end,end]*Ic


    # Get boundary condtion matrix
    M = P1*yR

    b = zeros(ComplexDF64, 3)
    b[3] = (2n+1)/r[end,end] 
    C2 = M \ b

    y = yR*C2

    # C2 = M[1:3,1:3] \ b[1:3]

    # y = yR[1:6,1:3]*C2


    k2 = y[5] - 1
    # println("Solid body k2 = ", k2)  

    k2_1[i] = k2

    Ediss1 = 21/2 * -imag(k2) * (ω0*R)^5/G * e^2
    Edot1[i] = Ediss1
end

# fig, axes = plt.subplots(ncols=8, figsize=(12,3.5))
fig, axes = plt.subplots(figsize=(4.5,3.5))
ytest = 0
disp = 0
vel = 0
for j in -7.3:1:-7.3
    for i in eachindex(ηs)
        #######################################################################
        η[3] = ηs[i] * (non_dim ? 1.0/(μ0*T) : 1.0)
        k[3] = 10.0^j * (non_dim ? 1.0/R^2 : 1.0)

        μc =  1im*ω*μ ./ (1im*ω .+ μ ./ η)

        # μc = μ

        
        # # Ic = get_Ic(r[end,1], ρ[1], g[end,1], μ[1], "liquid")

        # # Solve y!

        # B1 = get_B_product(r, ρ, g, μc, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k, 2, 2)[:, :, :, 2]
        # B2 = get_B_product(r, ρ, g, μc, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k, 3, 3)[:, :, :, 3]
        # B3 = get_B_product(r, ρ, g, μc, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k, 4, 4)[:, :, :, 4]
        # # display(B3*B2*B1 .- Bprod1[:,:,end,end])
        # # display()

        # Ic = get_Ic(r[end,1], ρ[1], g[end,1], μ[1], "liquid", 8, 4)

        # s1 = B1[:,:,end]*Ic

        # s1_start = copy(s1)
        # s1_start[7,4] = 1.0 # Set pore pressure to non-zero value
        # s2 = B2[:,:,end]*s1_start

        # s2_start = copy(s2)
        # s2_start[7:8, :] .= 0.0  # Set pore pressure and Darcy flux to zero value
        # s3 = B3[:,:,end]*s2
        ###

        

        global y = calculate_y(r, ρ, g, μc, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k)

        # println("Y test")
        # display(ytest)

        



        # # Construct boundary condition matrix
        # # M = zeros(ComplexDF64, 4,4)
        # M = zeros(precc, 4,4)

        # # Row 1 - Radial Stress
        # M[1, 1] = s3[3,1]
        # M[1, 2] = s3[3,2]
        # M[1, 3] = s3[3,3]
        # M[1, 4] = s3[3,4]

        # # Row 2 - Tangential Stress
        # M[2, 1] = s3[4,1]
        # M[2, 2] = s3[4,2]
        # M[2, 3] = s3[4,3]
        # M[2, 4] = s3[4,4]

        # # Row 3 - Potential Stress
        # M[3, 1] = s3[6,1]
        # M[3, 2] = s3[6,2]
        # M[3, 3] = s3[6,3]
        # M[3, 4] = s3[6,4]

        # #  Row 4 - Darcy flux (r = r_tp)
        # M[4, 1] = s2[8,1]
        # M[4, 2] = s2[8,2]
        # M[4, 3] = s2[8,3]
        # M[4, 4] = s2[8,4]

        # # b = zeros(ComplexF64, 4)
        # b = zeros(precc, 4)
        # b[3] = (2n+1)/r[end,end] 
        # C2 = M \ b

        # y = s3*C2

        # # display(y .- ytest[:,end,end])

        # # step1 = B1[:,:,end]*Ic*C2
        # # step1[7] = C2[4]
        # # step2 = B2[:,:,end]*step1
        # # step2[7:8] .= 0.0
        
        # # C3 = zeros(precc, 4)
        # # C3[1] = 1.0
        # # # println(size(yall), " ", size(r),size(B1))
        # # for k in 1:size(B1)[3]
        # #     y12 = B1[:,:,k]*Ic*C2
        # #     yall[:,k,2] .= y12[:] 
        # #     y23 = B2[:,:,k]*step1
        # #     yall[:,k,3] .= y23[:]
        # #     y34 = B3[:,:,k]*step2
        # #     yall[:,k,4] .= y34[:] 
        # #     # display(yall[:,k,4])
        # # end
        # # # display(yall[:,1,4])
        # # display(C2)

        yR = y[:,end,end]

        U22E =  7/8 * ω0^2*R^2*e * (non_dim ? ρ0/μ0 : 1)
        U22W = -1/8 * ω0^2*R^2*e * (non_dim ? ρ0/μ0 : 1)
        U20  = -3/2 * ω0^2*R^2*e * (non_dim ? ρ0/μ0 : 1)

        # global disp = get_displacement.(y[:,end,end-1], mag, 0.25π)
        sol_22  = get_solution(conj.(y), 2,  2, r, ρ, g, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k)
        sol_22c = get_solution(conj.(y), 2, -2, r, ρ, g, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k)
        sol_20  = get_solution(conj.(y), 2,  0, r, ρ, g, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k)

        global disp = U22E*sol_22[1] + U22W*sol_22c[1] + U20*sol_20[1] 
        # global disp = 0.5*(disp + conj.(disp))
        global vel = U22E*sol_22[2] + U22W*sol_22c[2] + U20*sol_20[2] 
        # global vel = 0.5*(vel + conj.(vel))


        k2 = yR[5] - 1

        k2_2[i] = k2    # println("Porous body k2 = ", k2)  
        println("k2 = ", yR[5] - 1)  
        println("h2 = ", -g[end]*yR[1] )
        Ediss2 = 21/2 * -imag(k2) * (ω0*R)^5/G * e^2
        Edot2[i] = Ediss2
        # println(Ediss1/1e9, " ", Ediss2/1e9)
        println(g[end,end])
        println(R^2 * ω0^2 * e / g[end,end] )

    end

    # ax2.loglog(ηs, Edot2/1e12, label="Melt, k=\$10^{$(j)}\$ m\$^2\$")
    # ax1.semilogx(ηs, real(k2_2))
end

# Now dimensionalise the solutions, if doing non-dimensional calculations
if non_dim
    disp[:] .*= R
    vel[:] .*= R / T
end

# Scale the solutions by (r/R)^2
for i in 1:size(r)[2]
    for j in 1:size(r)[1]-1
        scale = non_dim ? 1.0 : 1.0/R
        disp[:,:,:,j,i] .*= (r[j,i] * scale).^2
        vel[:,:,:,j,i] .*= (r[j,i] * scale ).^2
    end 
end

Eₗ_vol = zeros( (size(vel)[1], size(vel)[2], size(vel)[4]) )
Eₗ_area = zeros( (size(vel)[1], size(vel)[2]) )
println(size(vel))
display(k)
for i in 1:size(r)[1]-1
    # Dissipated energy per unit volume
    Eₗ_vol[:,:,i] = 0.5 * ηₗ[end-1]/k[end-1] * (abs.(vel[:,:,1,i,end-1]).^2 + abs.(vel[:,:,2,i,end-1]).^2 + abs.(vel[:,:,3,i,end-1]).^2)
    non_dim ? Eₗ_vol[:,:,i] *=  μ0 * T / R^2 : continue # dimensionalise ηₗ and k

    # Integrate across r to find dissipated energy per unit area
    dr = (r[i+1, end-1] - r[i, end-1]) * (non_dim ? R : 1.0)
    Eₗ_area[:,:] += Eₗ_vol[:,:,i] * dr
end
# To scale y functions for a given tida potential, use Eq. D6 in Rovira-Navarro et al., (2022)

# for i in 1:8
#     axes[i].plot(Double64.(real.(y[i,:,2:end][:])), R/1e3*Double64.(r[1:end-1, 2:end][:]), "-")
# end
# ax2.loglog(ηs, Edot1/1e12,"k--", label="No melt")
# println(size(disp))
# c = axes.contourf(disp[:,:,1,end,end])

# Eₗ_area = sum(Eₗ_vol*r[j,end-1], dims=3)[:,:,1]
# println(size(Eₗ_area))

c = axes.contourf(Eₗ_area)

PyPlot.colorbar(c)

# ax2.set_xlabel("Asthenosphere Solid Viscosity [Pa s]")
# ax2.set_ylabel("Tidal Heating Rate [TW]")

# ax2.axhspan((9.33-1.87)*1e13/1e12, (9.33+1.87)*1e13/1e12, alpha=0.5)

# ax1.semilogx(ηs, real(k2_1), "k--")

# # ax1.loglog(μs , k2_all, color="C1", label="Numerical")
# # ax1.loglog(μs, k2_an, "k:", label="Analytical (solid only)")
# # ax1.legend(frameon=false)
# ax1.set_xlabel("Asthenosphere Solid Viscosity [Pa s]")
# ax1.set_ylabel("k2 Tidal Love Number")

# ax2.legend(frameon=false, bbox_to_anchor=(1.0,0.5))

# ax2.set_xlim([ηs[1], ηs[end]])
# fig.subplots_adjust(wspace=.3)

# fig.savefig("io_porous.png", dpi=600, bbox_inches="tight")

show()


# # P3 = zeros(3, 8)
# # P3[1,3] = 1
# # P3[2,4] = 1
# # P3[3,6] = 1

# P4 = zeros(4, 8)
# P4[1,3] = 1
# P4[2,4] = 1
# P4[3,6] = 1
# P4[4,7] = 1

# Pl = zeros(Int64, 8,4)
# Pl[7,4] = 1

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
# # ybp = Bprod4[:,:,end,tp]*Ic + Bprod2[:,:,end,tp]*Pl

# # display(yR)
# # display(ytp)


# C1v = zeros(4)
# C2v = zeros(4)
# C3v = zeros(4)
# C4v = zeros(4)

# C1v[1] = 1
# C2v[2] = 1
# C3v[3] = 1
# C4v[4] = 1

# display(yR)
# display(yR *C1v)
# display(yR *C2v)
# display(yR *C3v)
# display(yR *C4v)

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
# b = zeros(ComplexF64, 4)
# b[3] = (2n+1)/r[end,end] 
# C2 = M \ b

# # C3 = M[1:3,1:3] \ b[1:3]

# y = yR*C2

# # display(ytp[:,1:3]*C3)
# # b[1] = (ytp[:,1:3]*C3)[3]
# # b[2] = (ytp[:,1:3]*C3)[4]
# # b[3] = (ytp[:,1:3]*C3)[6]


# # C4 = M \ b

# # # display(C2)
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

# k2 = y[5] - 1
# println("k2 = ", y[5] - 1)  
# println("h2 = ", -g[end]*y[1] )
# Ediss = 21/2 * -imag(k2) * (ω*r[end,end])^5/G * e^2
# println(Ediss/1e9)



