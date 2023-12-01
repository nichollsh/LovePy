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
    using DoubleFloats
    using AssociatedLegendrePolynomials
    include("SphericalHarmonics.jl")
    using .SphericalHarmonics

    export get_g, get_A!, get_A, get_B_product, get_Ic, get_B
    export expand_layers, set_G, calculate_y
    export get_displacement, get_darcy_velocity, get_solution

    prec = Double64 #BigFloat
    precc = ComplexDF64 #Complex{BigFloat}

    # prec = BigFloat
    # precc = Complex{BigFloat}

    G = prec(6.6743e-11)
    n = 2

    porous = false

    M = 6 + 2porous         # Matrix size: 6x6 if only solid material, 8x8 for two-phases
    nr = 150              # Number of sub-layers in each layer (TODO: change to an array)

    α = 0.95

    # Overwrite Gravitional constant for non-dimensional 
    # calculations
    function set_G(new_G)
        TidalLoveNumbers.G = new_G
    end


    function get_g(r, ρ)
        # g = zeros(Double64, size(r))
        # M = zeros(Double64, size(r))

        g = zeros(prec, size(r))
        M = zeros(prec, size(r))

        for i in 1:size(r)[2]
            M[2:end,i] = 4.0/3.0 * π .* diff(r[:,i].^3) .* ρ[i]
        end
        println(G)
        g[2:end,:] .= G*accumulate(+,M[2:end,:]) ./ r[2:end,:].^2
        g[1,2:end] = g[end,1:end-1]

        return g

    end

    function get_A(r, ρ, g, μ, κ)
        # A = zeros(ComplexDF64, 6, 6) 
        A = zeros(precc, 6, 6) 
        get_A!(A, r, ρ, g, μ, κ)
        return A
    end

    function get_A!(A::Matrix, r, ρ, g, μ, κ, λ=nothing)
        if isnothing(λ)
            λ = κ - 2μ/3
        end

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

        A[4,1] = -r_inv * (6κ*μ*r_inv*β_inv - ρ*g )
        A[4,2] = 2μ*r_inv^2 * (n*(n+1)*(1 + λ*β_inv) - 1.0 )
        A[4,3] = -r_inv * λ * β_inv # changed to match sabadini    
        A[4,4] = -3r_inv
        A[4,5] = ρ*r_inv

        A[5,1] = -4π * G * ρ
        A[5,5] = -(n+1)r_inv
        A[5,6] = 1.0

        A[6,1] = -4π*(n+1)*G*ρ*r_inv
        A[6,2] = -A[6,1]*n
        A[6,6] = (n-1)r_inv

    end

    # Method 2 for matrix propagator: two-phase flow
    function get_A(r, ρ, g, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k)
        # A = zeros(ComplexDF64, 8, 8)
        A = zeros(precc, 8, 8)
        get_A!(A, r, ρ, g, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k)
        return A
    end

    function get_A!(A::Matrix, r, ρ, g, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k)
        if α == 1
            κd = κ 
            κₛ = κd*1e5
            κu = κd + κₗ/ϕ
        else
            κₛ = κ
            κd = (1-α)κₛ
            if iszero(ϕ) || iszero(α)
                κu = 0.0
            else
                κu = κd + κₗ*κₛ*α^2 / (ϕ*κₛ + (α-ϕ)*κₗ)
            end

        end

        #κ is the bulk modulus of the solid! The drained bulk modulus
        # is (1-α)*κ


        mask = ϕ > 0.0

        mask = 1
        λ = κd .- 2μ/3
        κₛ = κ

        get_A!(A, r, ρ, g, μ, κd, λ) 

        r_inv = 1.0/r
        β_inv = 1.0/(2μ + λ)

        if !iszero(ϕ)
            A[1,7] = α * β_inv

            A[3,1] += 2g*ρₗ*α*r_inv * (-λ*β_inv +1) * mask           
            A[3,2] += n*(n+1)*g*ρₗ*α*r_inv * (λ*β_inv - 1) * mask      
            A[3,3] += β_inv*g*α*ρₗ * mask                             
            A[3,7] = (α*β_inv * (-4μ*r_inv + ρₗ*g*α) + g*ρₗ*(ϕ/κₗ + (α-ϕ)/κₛ)) * mask
            A[3,8] = 4π*G*ρₗ*ρ*k/(1im*ω*ηₗ) * mask
        
            A[4,7] = 2α*μ*r_inv * β_inv * mask
            
            A[5,8] = 4π*G*ρₗ*k / (1im * ω * ηₗ) * mask

            A[6,5] = -4π*n*(n+1)G*(ρₗ)^2*k*r_inv^2 / ( 1im*ω*ηₗ) * mask 
            A[6,7] = -4π*n*(n+1)G*ρₗ*k*r_inv^2 / ( 1im*ω*ηₗ) * mask
            A[6,8] = 4π*G*(n+1)*ρₗ*k*r_inv / (1im*ω*ηₗ) * mask
            
            A[7,1] = 4π*G*ρₗ*ρ * mask
            A[7,5] = ρₗ*(n+1)*r_inv * mask
            A[7,6] = -ρₗ * mask
            A[7,7] = -g*ρₗ / κₗ * mask
            A[7,8] = (1 - 4π*G*(ρₗ)^2*k / (1im*ω*ηₗ)) * mask
        
            A[8,1] = 2im*α*ω*ηₗ*r_inv/k * (1 - λ*β_inv) * mask
            A[8,2] = -1im*n*(n+1)*α*ω*ηₗ*r_inv/k * (1 - λ*β_inv) * mask
            A[8,3] = 1im*α*ω*ηₗ*β_inv / k * mask
            A[8,5] = n*(n+1)ρₗ * r_inv^2 * mask
            A[8,7] = ((n+1)n*r_inv^2 + 1im*ω*ηₗ/k * (α^2 *β_inv + ϕ /κₗ + (α-ϕ)/κₛ)) * mask
            A[8,8] = -2r_inv * mask
        end
        
    end

    function get_B_product(r, ρ, g, μ, κ)
       # Bprod = zeros(ComplexDF64, length(r), 8, 8)
    
    #    B = zeros(ComplexDF64, 6, 6)
       B = zeros(precc, 6, 6)
       B = I # Set B to the Identity matrix
       # Bprod[1,:,:] .= B[:,:]

       layer_num = size(r)[2]
       nr = size(r)[1]


    #    Bprod = zeros(ComplexDF64, (6, 6, nr-1, layer_num))
       Bprod = zeros(precc, (6, 6, nr-1, layer_num))

       for i in 2:layer_num # start at the top of the innermost layer
           r1 = r[1,i]
           for j in 1:nr-1
               r2 = r[j+1,i]
               g1 = g[j, i]
               g2 = g[j+1, i]
   
               B = get_B(r1, r2, g1, g2, ρ[i], μ[i], κ[i]) * B 

               Bprod[:,:,j,i] .= B[:,:]

               r1 = r2
           end
       end


        return Bprod
    end

    # second method: porous layer
    function get_B_product(r, ρ, g, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k, i1=2, iend=nothing)
        # Bprod = zeros(ComplexDF64, length(r), 8, 8)
    
        # B = zeros(ComplexDF64, 8, 8)
        B = zeros(precc, 8, 8)
        # B = I # Set B to the Identity matrix
        # B[7,7] = 0.0
        # Bprod[1,:,:] .= B[:,:]
        for i in 1:6
            B[i,i] = 1
        end

        # if starting from a porous layer, 
        # don't filter out y7 and y8 components
        if ϕ[i1]>0
            B[7,7] = 1
            B[8,8] = 1
        end

        layer_num = size(r)[2]
        nr = size(r)[1]


        # Bprod = zeros(ComplexDF64, (8, 8, nr-1, layer_num))
        Bprod = zeros(precc, (8, 8, nr-1, layer_num))

        r1 = r[1]
        if isnothing(iend)
            iend = layer_num
        end

        # Pδ = zeros(Int64, 8, 8)
        Pδ = zeros(BigInt, 8, 8)
        Pδ[1,1] = 1
        Pδ[2,2] = 1
        Pδ[3,3] = 1
        Pδ[4,4] = 1
        Pδ[5,5] = 1
        Pδ[6,6] = 1

        for i in i1:iend # start at the top of the innermost layer
            r1 = r[1,i]
            for j in 1:nr-1
                r2 = r[j+1,i]
                dr = r2 - r1
                rhalf = r1 + 0.5dr

                r2 = r[j+1, i]
                g1 = g[j, i]
                g2 = g[j+1, i]

                if ϕ[i] > 0 && j>1
                    Pδ[7,7] = 1
                    Pδ[8,8] = 1
                else
                    Pδ[7,7] = 0
                    Pδ[8,8] = 0
                end

                # In the first integration, don't filter out 
                # y7 and y8. Better way to do this? 
                if i==i1 && j==1
                    Pδ[7,7] = 1
                    Pδ[8,8] = 1
                end
    
                B = get_B(r1, r2, g1, g2, ρ[i], μ[i], κ[i], ω, ρₗ[i], κₗ[i], ηₗ[i], ϕ[i], k[i]) * B * Pδ
               
                Bprod[:,:,j,i] .= B[:,:]

                r1 = r2
            end
        end

        return Bprod
    end

    # Create R and S vectors?
    function set_sph_expansion(res=5.0)


    end

    function get_solution(y, n, m, r, ρ, g, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k, res=5.0)
        lons = deg2rad.(collect(0:res:360-0.001))'
        clats = deg2rad.(collect(0:res:180))
        cosTheta = cos.(clats)

        Y = m < 0 ? Ynmc(n,abs(m),clats,lons) : Ynm(n,abs(m),clats,lons)
        S = m < 0 ? Snmc(n,abs(m),clats,lons) : Snm(n,abs(m),clats,lons)
        
        # Y20 = Ynm(2,0,clats,lons) 
        # S20 = Snm(2,0,clats,lons)

        dY22dθ =
        dY20dθ = 

        dY22dϕ = 1im*m*Y
        dY20dϕ = 0.0

        d2Y22dθ2 = 
        d2Y20dθ2 =
        
        X22 = 
        X20 = 

        Z22 = 
        Z20 = 

        
        # Need to take conjugate in y here depending on whether the Fourier
        # transform is taken with exp(-iωt) or exp(iωt)
        # Also need negative sign due to the sign convention of the tidal potential
        # U22 = 0.5*mag*(Y22 * (7/8 * exp(-1im *ωt) ))
        # U22 = 0.5*mag*(Y22 * (7/8 * exp(-1im *ωt) - 1/8 *exp(1im *ωt) ))
        # U22 = 0.5*mag*(Y22 * (-1/8 *exp(1im *ωt) ))
        # U20 = 0.5*mag* -1.5(Y20 * exp(-1im * ωt) )
        # U20 .*= 0.0

        # U22_theta = 0.5*mag*(S22[1] * (7/8 * exp(-1im *ωt) - 1/8 *exp(1im *ωt) ))
        # U22_phi = 0.5*mag*(S22[2] * (7/8 * exp(-1im *ωt) - 1/8 *exp(1im *ωt) ))
        # U20_theta = 0.5*mag* -1.5(S20[1] * exp(-1im * ωt) )
        # U20_phi = 0.5*mag* -1.5(S20[2] * exp(-1im * ωt) )

        # mag_r = U22 + U20
        # mag_phi = U22_phi + U20_phi
        # mag_theta = U22_theta + U20_theta

        disp = zeros(ComplexF64, length(clats), length(lons), 3, size(r)[1]-1, size(r)[2])
        q_flux = zeros(ComplexF64, length(clats), length(lons), 3, size(r)[1]-1, size(r)[2])
        ϵrr = zero(disp)
        ϵθθ = zero(disp)
        ϵϕϕ = zero(disp)
        ϵrθ = zero(disp)
        ϵrϕ = zero(disp)
        ϵθϕ = zero(disp)
        ϵV = zero(disp)

        σrr = zero(disp)
        σθθ = zero(disp)
        σϕϕ = zero(disp)
        σrθ = zero(disp)
        σrϕ = zero(disp)
        σθϕ = zero(disp)
        # σV = zero(disp)
  
        for i in 2:size(r)[2] # Loop of layers
            ηₗr = ηₗ[i]
            ρₗr = ρₗ[i]
            ρr = ρ[i]
            kr  = k[i]
            κₗr = κₗ[i]
            κr = κ[i]
            μr = μ[i]
            ϕr = ϕ[i]
            for j in 1:size(r)[1]-1 # Loop over sublayers 
                # y[5,j,i], y[7,j,i], y[8,j,i], k[i], r[j,i], ηₗ[i], ρₗ[i]

                (y1, y2, y3, y4, y5, y6, y7, y8) = y[:,j,i]
                
                rr = r[j,i]
                gr = g[j,i]
                
                disp[:,:,:,j,i]   .= get_displacement(y1, y2, Y, S)
                if ϕ[i] > 0
                    q_flux[:,:,:,j,i] .= get_darcy_velocity(y5, y7, y8, kr, rr, ηₗr, ρₗr, Y, S)
                end

                A = get_A(rr, ρr, gr, μr, κr, ω, ρₗr, κₗr, ηₗr, ϕr, kr)
                dy1dr = A[1,:] * y[:,j,i]'
                dy2dr = A[2,:] * y[:,j,i]'

                # ϵrr = dy1dr * Y22
                # ϵθθ = 1/rr * ((y1 - 0.5n*(n+1)y2)Y22 + 0.5y2*X22)
                # ϵϕϕ = 1/rr * ((y1 - 0.5n*(n+1)y2)Y22 - 0.5y2*X22)
                # ϵrθ = 0.5 * (dy2dr + (y1 - y2)/rr) .* dY22dθ
                # ϵrϕ = 0.5 * (dy2dr + (y1 - y2)/rr) .* dY22dϕ .* 1.0 ./ sin.(clats) 
                # ϵθϕ = 0.5 * y2/rr * Z22
                # ϵV = dy1dr .+ 2/rr * y1 .- n*(n+1)/rr * y2

                # σrr = λ*ϵV * Y22 + 2μr*ϵrr - α*y7*Y22
                # σθθ = λ*ϵV * Y22 + 2μr*ϵθθ - α*y7*Y22
                # σϕϕ = λ*ϵV * Y22 + 2μr*ϵϕϕ - α*y7*Y22
                # σrθ = 2μr*ϵrθ
                # σrϕ = 2μr*ϵrϕ
                # σθϕ = 2μr*ϵθϕ
                
            end
        end

        
        

        return disp, q_flux#, (ϵrr, ϵθθ, ϵϕϕ, ϵrθ, ϵrϕ, ϵθϕ)
    end

    function get_displacement(y1, y2, Y, S)
        # displ_R =  mag_r * -conj.(y1)
        # displ_R .+= conj.(displ_R)

        # Conjugate y for eastward forcing, not westward
        displ_R =  Y * -y1
        # displ_R .+= conj.(displ_R)


        displ_theta = S[1] * -y2
        # displ_theta .+= conj.(displ_theta)

        displ_phi = S[2] * -y2
        # displ_phi .+= conj.(displ_phi)

        # Return drops the imaginary component, which should be zero anyway. Add check?
        # Radial component, theta component, phi component
        displ_vec = hcat(displ_R, displ_theta, displ_phi)  

        return reshape(displ_vec, (size(displ_R)[1], size(displ_R)[2], 3) )  
    end

    function get_darcy_velocity(y5, y7, y8, k, r, ηₗ, ρₗ, Y, S)
        q_R =  Y * y8 * -k/ηₗ
        # q_R .+= conj.(q_R)

        q_theta = S[1]  * -k/ηₗ * 1/r * ( y7 .+ ρₗ*y5) 
        # q_theta .+= conj.(q_theta)

        q_phi = S[2] * -k/ηₗ * 1/r * ( y7 .+ ρₗ*y5)
        # q_phi .+= conj.(q_phi)

        # Return drops the imaginary component, which should be zero anyway. Add check?
        # Radial component, theta component, phi component

        q_vec = hcat(q_R, q_theta, q_phi)  

        return reshape(q_vec, (size(q_R)[1], size(q_R)[2], 3) )   
    end


    # function get_displacement(y1, y2)
    #     # lons = deg2rad.(collect(0:res:360-0.001))'
    #     # clats = deg2rad.(collect(0:res:180))
    #     # cosTheta = cos.(clats)

    #     # Y22 = Ynm(2,2,clats,lons)
    #     # S22 = Snm(2,2,clats,lons)
        
    #     # Y20 = Ynm(2,0,clats,lons) 
    #     # S20 = Snm(2,0,clats,lons)
        
    #     # # Need to take conjugate in y here depending on whether the Fourier
    #     # # transform is taken with exp(-iωt) or exp(iωt)
    #     # # Also need negative sign due to the sign convention of the tidal potential
    #     # # U22 = -0.5*ω^2*R^2*e*(Y22 * (7/8 * exp(-1im *ωt) )) * conj.(y[1,end,end])
    #     # U22 = 0.5*mag*(Y22 * (7/8 * exp(-1im *ωt) - 1/8 *exp(1im *ωt) ))
    #     # U20 = 0.5*mag* -1.5(Y20 * exp(-1im * ωt) )

    #     # U22_theta = 0.5*mag*(S22[1] * (7/8 * exp(-1im *ωt) - 1/8 *exp(1im *ωt) ))
    #     # U22_phi = 0.5*mag*(S22[2] * (7/8 * exp(-1im *ωt) - 1/8 *exp(1im *ωt) ))
    #     # U20_theta = 0.5*mag* -1.5(S20[1] * exp(-1im * ωt) )
    #     # U20_phi = 0.5*mag* -1.5(S20[2] * exp(-1im * ωt) )

    #     displ_R22 =  U22 * -conj.(y[1])
    #     displ_R22 .+= conj.(displ_R22)

    #     displ_R20 = U20 * -conj.(y[1])
    #     displ_R20 .+= conj.(displ_R20) 

    #     displ_S22_theta = U22_theta * -conj.(y[2])
    #     displ_S22_theta .+= conj.(displ_S22_theta)

    #     displ_S22_phi = U22_phi * -conj.(y[2])
    #     displ_S22_phi .+= conj.(displ_S22_phi)

    #     displ_S20_theta = U20_theta * -conj.(y[2])
    #     displ_S20_theta .+= conj.(displ_S20_theta)

    #     displ_S20_phi = U20_phi * -conj.(y[2])
    #     displ_S20_phi .+= conj.(displ_S20_phi)

    #     # Return drops the imaginary component, which should be zero anyway. Add check?
    #     # Radial component, theta component, phi component
    #     displ_vec = hcat(real(displ_R22 + displ_R20), real(displ_S22_theta + displ_S20_theta), real(displ_S22_phi + displ_S20_phi))  
        
    #     return reshape(displ_vec, (length(clats), length(lons), 3) )  
    # end

    function get_strain(y)
    end

    function get_stress(y)
    end

    # function get_darcy_velocity(y, mag, k, r, ηₗ, ρₗ, ωt=0.0, res=5.0,n=2, m=2)
    #     lons = deg2rad.(collect(0:res:360-0.001))'
    #     clats = deg2rad.(collect(0:res:180))
    #     cosTheta = cos.(clats)

    #     Y22 = Ynm(2,2,clats,lons)
    #     S22 = Snm(2,2,clats,lons)
        
    #     Y20 = Ynm(2,0,clats,lons) 
    #     S20 = Snm(2,0,clats,lons)
        
    #     # Need to take conjugate in y here depending on whether the Fourier
    #     # transform is taken with exp(-iωt) or exp(iωt)
    #     # Also need negative sign due to the sign convention of the tidal potential
    #     # U22 = -0.5*ω^2*R^2*e*(Y22 * (7/8 * exp(-1im *ωt) )) * conj.(y[1,end,end])
    #     U22 = -0.5*mag*(Y22 * (7/8 * exp(-1im *ωt) - 1/8 *exp(1im *ωt) ))
    #     U20 = -0.5*mag* -1.5(Y20 * exp(-1im * ωt) )

    #     U22_theta = -0.5*mag*(S22[1] * (7/8 * exp(-1im *ωt) - 1/8 *exp(1im *ωt) ))
    #     U22_phi = -0.5*mag*(S22[2] * (7/8 * exp(-1im *ωt) - 1/8 *exp(1im *ωt) ))
    #     U20_theta = -0.5*mag* -1.5(S20[1] * exp(-1im * ωt) )
    #     U20_phi = -0.5*mag* -1.5(S20[2] * exp(-1im * ωt) )

    #     q_R22 =  U22 * conj.(y[8,end-10,end-1]) * -k[end-1]/ηₗ[end-1]
    #     q_R22 .+= conj.(q_R22)

    #     q_R20 = U20 * conj.(y[8,end-10,end-1]) * -k[end-1]/ηₗ[end-1]
    #     q_R20 .+= conj.(q_R20) 

    #     q_S22_theta = U22_theta  * -k[end-1]/ηₗ[end-1] * 1/r * ( conj.(y[7,end,end-1]) .+ ρₗ[end-1]*conj.(y[5,end,end-1])) 
    #     q_S22_theta .+= conj.(q_S22_theta)

    #     q_S22_phi = U22_phi * -k[end-1]/ηₗ[end-1] * 1/r * ( conj.(y[7,end,end-1]) .+ ρₗ[end-1]*conj.(y[5,end,end-1]))
    #     q_S22_phi .+= conj.(q_S22_phi)

    #     q_S20_theta = U20_theta * -k[end-1]/ηₗ[end-1] * 1/r * ( conj.(y[7,end,end-1]) .+ ρₗ[end-1]*conj.(y[5,end,end-1]))
    #     q_S20_theta .+= conj.(q_S20_theta)

    #     q_S20_phi = U20_phi * -k[end-1]/ηₗ[end-1] * 1/r * ( conj.(y[7,end,end-1]) .+ ρₗ[end-1]*conj.(y[5,end,end-1]))
    #     q_S20_phi .+= conj.(q_S20_phi)

    #     # Return drops the imaginary component, which should be zero anyway. Add check?
    #     # Radial component, theta component, phi component

    #     q_vec = hcat(real.(q_R22 + q_R20), real(q_S22_theta + q_S20_theta), real(q_S22_phi + q_S20_phi))  

    #     return reshape(q_vec, (length(clats), length(lons), 3) )   
    # end

    

    function calculate_y(r, ρ, g, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k, core="liquid")
        porous_layer = ϕ .> 0.0

        sum(porous_layer) > 1.0 && error("Can only handle one porous layer for now!")

        nlayers = size(r)[2]
        nsublayers = size(r)[1]

        Ic = get_Ic(r[end,1], ρ[1], g[end,1], μ[1], core, 8, 4)

        y1_4 = zeros(precc, 8, 4, size(r)[1]-1, size(r)[2]) # Four linearly independent y solutions
        y = zeros(precc, 8, size(r)[1]-1, size(r)[2])

        for i in 2:nlayers
            B = get_B_product(r, ρ, g, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k, i, i)[:, :, :, i]

            # Determine starting y vector for layer i
            if i==2
                y_start = copy(Ic)
            else
                y_start = y1_4[:,:,end,i-1]
            end

            # Modify starting vector for porous layer
            if porous_layer[i] && !porous_layer[i-1]
                y_start[7,4] = 1.0
            elseif !porous_layer[i]
                y_start[7:8, :] .= 0.0
            end

            for j in 1:nsublayers-1
                y1_4[:,:,j,i] = B[:,:,j] * y_start #y1_4[:,:,1,i]
            end

        end

        M = zeros(precc, 4,4)

        # Row 1 - Radial Stress
        M[1, :] .= y1_4[3,:,end,end]

        # Row 2 - Tangential Stress
        M[2, :] .= y1_4[4,:,end,end]
    
        # Row 3 - Potential Stress
        M[3, :] .= y1_4[6,:,end,end]
        
        #  Row 4 - Darcy flux (r = r_tp)
        for i in 2:nlayers
            if porous_layer[i]
                M[4, :] .= y1_4[8,:,end,i]
            end
        end
        
        b = zeros(precc, 4)
        b[3] = (2n+1)/r[end,end] 
        C = M \ b

        # y = y1_4[:,:,end,end]*C

        for i in 2:nlayers
            for j in 1:nsublayers-1
                y[:,j,i] = y1_4[:,:,j,i]*C
            end
        end

        # step1 = B1[:,:,end]*Ic*C2
        # step1[7] = C2[4]
        # step2 = B2[:,:,end]*step1
        # step2[7:8] .= 0.0
        
        # C3 = zeros(precc, 4)
        # C3[1] = 1.0
        # # println(size(yall), " ", size(r),size(B1))
        # for k in 1:size(B1)[3]
        #     y12 = B1[:,:,k]*Ic*C2
        #     yall[:,k,2] .= y12[:] 
        #     y23 = B2[:,:,k]*step1
        #     yall[:,k,3] .= y23[:]
        #     y34 = B3[:,:,k]*step2
        #     yall[:,k,4] .= y34[:] 
        #     # display(yall[:,k,4])
        # end

        # B1 = get_B_product(r, ρ, g, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k, 2, 2)[:, :, :, 2]
        # B2 = get_B_product(r, ρ, g, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k, 3, 3)[:, :, :, 3]
        # B3 = get_B_product(r, ρ, g, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k, 4, 4)[:, :, :, 4]

        

        return y

        # s1_start = copy(s1)
        # s1_start[7,4] = 1.0 # Set pore pressure to non-zero value
        # s2 = B2[:,:,end]*s1_start

        # s2_start = copy(s2)
        # s2_start[7:8, :] .= 0.0  # Set pore pressure and Darcy flux to zero value
        # s3 = B3[:,:,end]*s2

    end

    function get_B(r1, r2, g1, g2, ρ, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k)
        dr = r2 - r1
        rhalf = r1 + 0.5dr
        
        ghalf = g1 + 0.5*(g2 - g1)

        A1 = get_A(r1, ρ, g1, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k)
        Ahalf = get_A(rhalf, ρ, ghalf, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k)
        A2 = get_A(r2, ρ, g2, μ, κ, ω, ρₗ, κₗ, ηₗ, ϕ, k)
        
        k1 = dr * A1 
        k2 = dr * Ahalf * (I + 0.5k1)
        k3 = dr * Ahalf * (I + 0.5k2)
        k4 = dr * A2 * (I + k3) 

        B = (I + 1.0/6.0 * (k1 + 2k2 + 2k3 + k4))

        return B
    end

    function get_B(r1, r2, g1, g2, ρ, μ, κ)
        dr = r2 - r1
        rhalf = r1 + 0.5dr
        
        ghalf = g1 + 0.5*(g2 - g1)

        A1 = get_A(r1, ρ, g1, μ, κ)
        Ahalf = get_A(rhalf, ρ, ghalf, μ, κ)
        A2 = get_A(r2, ρ, g2, μ, κ)
        
        k1 = dr * A1 
        k2 = dr * Ahalf * (I + 0.5k1)
        k3 = dr * Ahalf * (I + 0.5k2)
        k4 = dr * A2 * (I + k3) 

        B = (I + 1.0/6.0 * (k1 + 2k2 + 2k3 + k4))

        return B
    end

    function get_Ic(r, ρ, g, μ, type, M=6, N=3)
        # Ic = zeros(Double64, M, N)
        Ic = zeros(prec, M, N)

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
        # rs = zeros(Double64, (nr+1, length(r)-1))
        rs = zeros(prec, (nr+1, length(r)-1))
        
        for i in 1:length(r)-1
            rfine = LinRange(r[i], r[i+1], nr+1)
            rs[:, i] .= rfine[1:end] 
        end
    
        return rs
    end

    
end