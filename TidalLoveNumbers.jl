# Julia code to calculate tidal deformation, Love numbers, and heating
# Author: H. Hay

# μ: solid shear modulus
# ρ: solid density 
# ρₗ: liquid density 
# κ: solid bulk modulus 
# κₗ: liquid bulk modulus 
# α: Biot's constant 
# λ: Lame's First Parameter
# η: solid viscosity 
# ηₗ: liquid viscosity 
# g: gravity 
# ϕ: porosity 
# k: permeability 
# ω: rotation rate


module TidalLoveNumbers

    using LinearAlgebra

    export get_g, get_A!, get_A, get_B_product, get_Ic
    export expand_layers

    G = 6.67408e-11
    n = 2

    porous = false

    M = 6 + 2porous         # Matrix size: 6x6 if only solid material, 8x8 for two-phases
    nr = 20               # Number of sub-layers in each layer (TODO: change to an array)

    α = 0.1             
    k = 1e-10
    # ω = 2*2.05e-5
    # ϕ = 0.1
    # ηₗ = 1e-2
    # η = 1e21
    ρₗ = 1e3

    function get_g(r, ρ)
        g = zero(r)
        M = zero(r)

        for i in 1:size(r)[2]
            M[2:end,i] = 4.0/3.0 * π .* diff(r[:,i].^3) .* ρ[i]
        end

        g[2:end,:] .= G*accumulate(+,M[2:end,:]) ./ r[2:end,:].^2
        g[1,2:end] = g[end,1:end-1]

        return g

    end

    function get_A(r, ρ, g, μ, κ, η, ω)
        A = zeros(ComplexF64, 6, 6) 
        get_A!(A, r, ρ, g, μ, κ, η, ω)
        return A
    end

    function get_A!(A::Matrix, r, ρ, g, μ, κ, η, ω)
        μ =  1im*ω*μ / (1im*ω + μ/η)    # Complex shear modulus for a Maxwell rheology

        λ = κ .- 2μ/3

        r_inv = 1.0/r
        β_inv = 1.0/(2μ + λ)

        A[1,1] = -2λ * β_inv * r_inv
        A[1,2] = n*(n+1) * λ * β_inv * r_inv
        A[1,3] = β_inv

        A[2,1] = -r_inv
        A[2,2] = r_inv
        A[2,4] = 1.0 / μ

        A[3,1] = 4r_inv * (3κ*μ*r_inv*β_inv - ρ*g) #+ 
                #  porous*(2g*ρₗ*α*r_inv * (-λ*β_inv +1))             # two-phase
        A[3,2] = -n*(n+1)*r_inv * (6κ*μ*r_inv*β_inv - ρ*g ) #+
                #  porous*(n*(n+1)*g*ρₗ*α*r_inv * (λ*β_inv - 1))      # two-phase
        A[3,3] = β_inv * (-4μ*r_inv )
                #  + porous*g *α*ρₗ)                                  # two-phase                                                
        A[3,4] = n*(n+1)*r_inv
        A[3,5] = -ρ * (n+1)*r_inv
        A[3,6] = ρ
        # if porous
        #     A[3,7] = α*β_inv * 
        #              (-4μ*r_inv + ρₗ*g*α) + g*ρₗ*(ϕ/κₗ + (α-ϕ)/κ)
        #     A[3,8] = 4π*G*ρₗ*ρ*k/(1im*ω*ηₗ)
        # end

        A[4,1] = -r_inv * (6κ*μ*r_inv*β_inv - ρ*g )
        A[4,2] = 2μ*r_inv^2 * (n*(n+1)*(1 + λ*β_inv) - 1.0 )
        A[4,3] = -r_inv * λ * β_inv # changed to match sabadini    
        A[4,4] = -3r_inv
        A[4,5] = ρ*r_inv
        # if porous
        #     A[4,7] = 2α*μ*r_inv / (2μ + λ)
        # end

        A[5,1] = -4π * G * ρ
        A[5,5] = -(n+1)r_inv
        A[5,6] = 1.0
        # if porous
        #     A[5,8] = 4π*G*ρₗ*k / (1im * ω * ηₗ)
        # end

        A[6,1] = -4π*(n+1)*G*ρ*r_inv
        A[6,2] = -A[6,1]*n
        # A[6,5] = porous * (-4π*n*(n+1)G*(ρₗ)^2*k*r_inv^2 / ( 1im*ω*ηₗ)) # two-phase
        A[6,6] = (n-1)r_inv
        # if porous                                                       # two-phase
        #     A[6,7] = -4π*n(n+1)G*ρₗ*k*r_inv^2 / ( 1im*ω*ηₗ)
        #     A[6,8] = 4π*G*(n+1)*ρₗ*k*r_inv / (1im*ω*ηₗ)
        # end

        # if porous                                                       # two-phase
        #     A[7,1] = 4π*G*ρₗ*ρ
        #     A[7,5] = ρₗ*(n+1)*r_inv
        #     A[7,6] = ρₗ
        #     A[7,7] = -g*ρₗ / κₗ
        #     A[7,8] = 1 - 4π*G*(ρₗ)^2*k / (1im*ω*ηₗ)
        
        #     A[8,1] = 2im*α*ω*ηₗ*r_inv/k * (1 - λ*β_inv)
        #     A[8,2] = -1im*n*(n+1)*α*ω*ηₗ*r_inv/k * (1 - λ*β_inv)
        #     A[8,3] = 1im*α*ω*ηₗ*β_inv / κ 
        #     A[8,5] = n*(n+1)ρₗ * r_inv^2
        #     A[8,7] = (n+1)n*r_inv^2 + 1im*ω*ηₗ/k * 
        #              (α^2 *β_inv + ϕ /κₗ + (α-ϕ)/κ)
        #     A[8,8] = -2r_inv
        # end

        # return A
        # display(A)
    end

    # Method 2 for matrix propagator: two-phase flow
    function get_A(r, ρ, g, μ, κ, ρₗ, α, κₗ, ω, ηₗ, ϕ)
        A = zeros(ComplexF64, 8, 8)
        get_A!(A, r, ρ, g, μ, κ, ρₗ, α, κₗ, ω, ηₗ, ϕ)
        return A
    end

    function get_A!(A::Matrix, r, ρ, g, μ, κ, ρₗ, α, κₗ, ω, ηₗ, ϕ)
        get_A!(A[1:6, 1:6], r, ρ, g, μ, κ, η, ω)  

        λ = κ .- 2μ/3

        r_inv = 1.0/r
        β_inv = 1.0/(2μ + λ)

        mask = ϕ > 0.0

        A[3,1] += 2g*ρₗ*α*r_inv * (-λ*β_inv +1) * mask           
        A[3,2] += n*(n+1)*g*ρₗ*α*r_inv * (λ*β_inv - 1) * mask      
        A[3,3] += β_inv*g*α*ρₗ * mask                             
        A[3,7] = (α*β_inv * (-4μ*r_inv + ρₗ*g*α) + g*ρₗ*(ϕ/κₗ + (α-ϕ)/κ)) * mask
        A[3,8] = 4π*G*ρₗ*ρ*k/(1im*ω*ηₗ) * mask
    
        A[4,7] = 2α*μ*r_inv / (2μ + λ) * mask
        
        A[5,8] = 4π*G*ρₗ*k / (1im * ω * ηₗ) * mask

        A[6,5] = -4π*n*(n+1)G*(ρₗ)^2*k*r_inv^2 / ( 1im*ω*ηₗ) * mask 
        A[6,7] = -4π*n*(n+1)G*ρₗ*k*r_inv^2 / ( 1im*ω*ηₗ) * mask
        A[6,8] = 4π*G*(n+1)*ρₗ*k*r_inv / (1im*ω*ηₗ) * mask
        
        A[7,1] = 4π*G*ρₗ*ρ * mask
        A[7,5] = ρₗ*(n+1)*r_inv * mask
        A[7,6] = ρₗ * mask
        A[7,7] = -g*ρₗ / κₗ * mask
        A[7,8] = 1 - 4π*G*(ρₗ)^2*k / (1im*ω*ηₗ) * mask
    
        A[8,1] = 2im*α*ω*ηₗ*r_inv/k * (1 - λ*β_inv) * mask
        A[8,2] = -1im*n*(n+1)*α*ω*ηₗ*r_inv/k * (1 - λ*β_inv) * mask
        A[8,3] = 1im*α*ω*ηₗ*β_inv / κ * mask
        A[8,5] = n*(n+1)ρₗ * r_inv^2 * mask
        A[8,7] = ((n+1)n*r_inv^2 + 1im*ω*ηₗ/k * (α^2 *β_inv + ϕ /κₗ + (α-ϕ)/κ)) * mask
        A[8,8] = -2r_inv * mask

        return A
        
    end

    function get_B_product(r, ρ, g, μ, κ, η, ω)
        
    
        B = zeros(ComplexF64, 6, 6)
        B = I # Set B to the Identity matrix
        # Bprod[1,:,:] .= B[:,:]

        #
        layer_num = size(r)[2]
        nr = size(r)[1]

        Bprod = zeros(ComplexF64, (6, 6, nr-1, layer_num))

        for i in 2:layer_num # start at the top of the innermost layer
            r1 = r[1,i]
            for j in 1:nr-1
                r2 = r[j+1,i]
                dr = r2 - r1
                rhalf = r1 + 0.5dr
                
                ghalf = g[j, i] + 0.5*(g[j+1, i] - g[j, i])

                A1 = get_A(r1, ρ[i], g[j,i], μ[i], κ[i], η[i], ω)
                Ahalf = get_A(rhalf, ρ[i], ghalf, μ[i], κ[i], η[i], ω)
                A2 = get_A(r2, ρ[i], g[j+1,i], μ[i], κ[i], η[i], ω)
                
                
                k1 = dr * A1 
                k2 = dr * Ahalf * (I + 0.5k1)
                k3 = dr * Ahalf * (I + 0.5k2)
                k4 = dr * A2 * (I + k3) 

                B = (I + 1.0/6.0 * (k1 + 2k2 + 2k3 + k4)) * B

                # display(B)

                Bprod[:,:,j,i] .= B[:,:]

                r1 = r2
            end
        end


        # for i in eachindex(r[1:end-1])
            
        #     r2 = r[i+1]
        #     dr = r2 - r1
        #     rhalf = r1 + 0.5dr
            
        #     ghalf = g[i] + 0.5*(g[i+1] - g[i])

        #     A1 = get_A(r1, ρ[i+1], g[i], μ[i+1], κ[i+1], η[i+1], ω)
        #     Ahalf = get_A(rhalf, ρ[i+1], ghalf, μ[i+1], κ[i+1], η[i+1], ω)
        #     A2 = get_A(r2, ρ[i+1], g[i+1], μ[i+1], κ[i+1], η[i+1], ω)
            
        #     k1 = dr * A1 
        #     k2 = dr * Ahalf * (I + 0.5k1)
        #     k3 = dr * Ahalf * (I + 0.5k2)
        #     k4 = dr * A2 * (I + k3) 

        #     B = (I + 1.0/6.0 * (k1 + 2k2 + 2k3 + k4)) * B

        #     Bprod[i+1,:,:] .= B[:,:]

        #     r1 = r2
        # end

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
        else # incompressible solid core
            # First column
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

    function expand_layers(r)#, ρ, μ, κ, η)
        # rs = zeros(Float64, (length(r)-1)*nr - length(r) + 2)
        rs = zeros(Float64, (nr+1, length(r)-1))
        # display(rs)
        # display(rs)


        # ρs = zero(rs)
        # μs = zeros(ComplexF64, size(rs))
        # κs = zero(rs)
        # ηs = zero(rs)

        

        # rindex = 1
        # ρs[1] = ρ[1]
        # μs[1] = μ[1]
        # κs[1] = κ[1]
        # ηs[1] = η[1]
        
            for i in 1:length(r)-1
                rfine = LinRange(r[i], r[i+1], nr+1)
                rs[:, i] .= rfine[1:end] 
                # ρs[:,i] .= ρ[i]
                # μs[rindex+1:rindex + length(rfine)  ] .= μ[i+1]
                # κs[rindex+1:rindex + length(rfine)  ] .= κ[i+1]
                # ηs[rindex+1:rindex + length(rfine)  ] .= η[i+1]
                
                # rindex = rindex + length(rfine)
            end
        
        # rs[end] = r[end]
        # display(rs)
        # display(ρs)
        
        return rs#, ρs, μs, κs, ηs
    end

    
end