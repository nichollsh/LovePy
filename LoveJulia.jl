# Julia code to calculate tidal deformation, Love numbers, and heating
# Author: H. Hay

using LinearAlgebra

G = 6.67408e-11
n = 2

porous = false

M = 6 + 2porous
nr = 60

function get_g(r, ρ)
    g = zeros(Float64, length(r))
    M = zeros(Float64, length(r))

    M[2:end] = 4.0/3.0 * π .* diff(r.^3) .* ρ[2:end]
    M[1] = 4.0/3.0 * π * r[1]^3 * ρ[1]

    g .= G*cumsum(M) ./ r.^2

    return g

end

function get_A(r, ρ, g, μ, κ)
    A = zeros(Float64, 6, 6)

    λ = κ .- 2μ/3

    r_recip = 1.0/r
    β_recip = 1.0/(2μ + λ)

    A[1,1] = -2λ * β_recip * r_recip
    A[1,2] = n*(n+1) * λ * β_recip * r_recip
    A[1,3] = β_recip

    A[2,1] = -r_recip
    A[2,2] = r_recip
    A[2,4] = 1.0 / μ

    A[3,1] = 4r_recip * (3κ*μ*r_recip*β_recip - ρ*g)
    A[3,2] = -n*(n+1)*r_recip * (6κ*μ*r_recip*β_recip - ρ*g )
    A[3,3] = β_recip * (-4μ*r_recip )
    A[3,4] = n*(n+1)*r_recip
    A[3,5] = -ρ * (n+1)*r_recip
    A[3,6] = ρ

    A[4,1] = -r_recip * (6κ*μ*r_recip*β_recip - ρ*g )
    A[4,2] = 2μ*r_recip^2 * (n*(n+1)*(1 + λ * β_recip) - 1.0 )
    A[4,3] = -r_recip * λ * β_recip # changed to match sabadini    
    A[4,4] = -3r_recip
    A[4,5] = ρ*r_recip

    A[5,1] = -4π * G * ρ
    A[5,5] = -(n+1)r_recip
    A[5,6] = 1.0

    A[6,1] = -4π*(n+1)*G*ρ*r_recip
    A[6,2] = -A[6,1]*n
    A[6,6] = (n-1)r_recip

   return A
end

function get_B_product(r, ρ, g, μ, κ)
    Bprod = zeros(length(r), M, M)
   
    B = zeros(M, M)
    B = I # Set B to the Identity matrix
    # Bprod[1,:,:] .= B[:,:]

    r1 = r[1]
    for i in eachindex(r[1:end-1])
        
        r2 = r[i+1]
        dr = r2 - r1
        rhalf = r1 + 0.5dr
        
        ghalf = g[i] + 0.5*(g[i+1] - g[i])

        A1 = get_A(r1, ρ[i+1], g[i], μ[i+1], κ[i+1])
        Ahalf = get_A(rhalf, ρ[i+1], ghalf, μ[i+1], κ[i+1])
        A2 = get_A(r2, ρ[i+1], g[i+1], μ[i+1], κ[i+1])
        
        k1 = dr * A1 
        k2 = dr * Ahalf * (I + 0.5k1)
        k3 = dr * Ahalf * (I + 0.5k2)
        k4 = dr * A2 * (I + k3) 

        B = (I + 1.0/6.0 * (k1 + 2k2 + 2k3 + k4)) * B

        Bprod[i+1,:,:] .= B[:,:]

        r1 = r2
    end

    return Bprod
end

function get_Ic(r, ρ, g, μ, type)
    Ic = zeros(Float64, M, 3)

    if type=="liquid"
        Ic[1,1] = -r^n / g
        Ic[1,3] = 1.0
        Ic[2,2] = 1.0
        Ic[3,3] = g*ρ
        Ic[5,1] = r^n
        Ic[6,1] = 2(n-1)*r^(n-1)
        Ic[6,3] = 4π * G * ρ 
    else
        Ic[1, 1] = n*r^( n+1 ) / ( 2*( 2n + 3) )
        Ic[2, 1] = ( n+3 )*r^( n+1 ) / ( 2*( 2n+3 ) * ( n+1 ) )
        Ic[3, 1] = ( n*ρ*g*r + 2*( n^2 - n - 3)*μ ) * r^n / ( 2*( 2n + 3) )
        Ic[4, 1] = n *( n+2 ) * μ * r^n / ( ( 2n + 3 )*( n+1 ) )
        Ic[6, 1] = 2π*G*ρ*n*r^( n+1 ) / ( 2n + 3 )

        # Second column

        Ic[1, 2] = r^( n-1 )
        Ic[2, 2] = r^( n-1 ) / n
        Ic[3, 2] = ( ρ*g*r + 2*( n-1 )*μ ) * r^( n-2 )
        Ic[4, 2] = 2*( n-1 ) * μ * r^( n-2 ) / n
        Ic[6, 2] = 4π*G*ρ*r^( n-1 )

        # Third column

        Ic[3, 3] = -ρ * r^n
        Ic[5, 3] = -r^n
        Ic[6, 3] = -( 2n + 1) * r^( n-1 )

    end

    return Ic
end

function expand_layers(r, ρ, μ, κ)
    rs = zeros(Float64, (length(r)-1)*nr - length(r) + 2)
    ρs = zero(rs)
    μs = zero(rs)
    κs = zero(rs)

    rindex = 1
    ρs[1] = ρ[1]
    μs[1] = μ[1]
    κs[1] = κ[1]
    for i in 1:length(r)-1
        rfine = LinRange(r[i], r[i+1], nr)[1:end-1]
        rs[rindex:rindex + length(rfine) - 1] .= rfine[1:end] 
        ρs[rindex+1:rindex + length(rfine)  ] .= ρ[i+1]
        μs[rindex+1:rindex + length(rfine)  ] .= μ[i+1]
        κs[rindex+1:rindex + length(rfine)  ] .= κ[i+1]
        
        rindex = rindex + length(rfine)
    end
    rs[end] = r[end]
    
    return rs, ρs, μs, κs
end

#enceladus test model:

ρ = [5376.81502017327, 3300, 1000, 1000, 940]
r = [720, 1441, 1541, 1556, 1561] .* 1e3
μ = [40, 40, 1.0, 3.5, 3.5] .* 1e9
κ = [5000e9, 5000e9, 5000e9, 5000e9, 5000e9]

# ρ = [7640, 7640, 3300, 3300, 3300, 3000]
# r = [20e3, 700e3, 1500e3, 1781.6e3, 1801.6e3, 1821.6e3]
# μ = [60e9, 60e9, 60e9, 60e9, 1e1, 60e9]
# κ = [5000e9, 5000e9, 5000e9, 5000e9, 5000e9, 5000e9]



# ρ = [5000, 3000, 3000]
# r = [130e3, 200e3, 2500e3]
# μ = [40e9, 40e9, 40e9]
# κ = [5000e9, 5000e9, 5000e9]

r, ρ, μ, κ = expand_layers(r, ρ, μ, κ)

# display(ρ)
# display(r)

# r = collect(1e3:1e3:2500e3)
# N = length(r)
# ρ = ones(N) .* 3000

# μ = ones(N) .* 3e9
# κ = ones(N) .* 5000000e9

# display(hcat(r, ρ, μ, κ))


g = get_g(r, ρ)

Ic = get_Ic(r[1], ρ[1], g[1], μ[1], "solid")

# display(Ic)

Bprod = get_B_product(r, ρ, g, μ, κ)[end,:,:]
y1 = Bprod*Ic
# display(y1)

P1 = zeros(Float64, 3, 6)
P1[1,3] = 1
P1[2,4] = 1
P1[3,6] = 1

b = zeros(Float64, 3)
b[3] = (2n+1)/r[end] 

C = (P1 * y1) \ b

# Seems to be equivalent to using P1 and P2 
y = C[1]*y1[:,1] + C[2]*y1[:,2] + C[3]*y1[:,3]

# display(y)

P2 = zeros(Float64, 3, 6)
P2[1,1] = 1
P2[2,2] = 1
P2[3,5] = 1


y1 = P1 * Bprod * Ic * C 
y2 = P2 * Bprod * Ic * C 


# display(y1)
# display(y2)
# display(b)


println( y[5] - 1,  ", ",3/2 / (1 + 19/2 * μ[end] / (ρ[end]*g[end]*r[end])  ))
println( -g[end]*y[1], ", ", 5/2 / (1 + 19/2 * μ[end] / (ρ[end]*g[end]*r[end])  ))
println( -g[end]*y[2]  )

