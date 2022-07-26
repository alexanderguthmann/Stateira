function radial_func(r)
    a0 = 5.29177210903E-11
    rs = 1.0E-10
    Eh = 4.3597447222071E-18
    h = 6.62607015E-34
    c0 = 2.99792458E08
    R_SR_S1 = 0.522E-09
    R_LR_S1 = 1.200E-09
    A_SR_S1 = -0.500680370E05
    B_SR_S1 = 4.34413885E11
    N_S1 = 7
    --coefficients for long range potential
    C6_S1 = 5693.7056 * Eh * a0^6
    C8_S1 = 796487.36 * Eh * a0^8
    C10_S1 = 95332817 * Eh * a0^10
    Aex_S1 = 0.37664685E05
    gamma_S1 = 5.427916
    beta_S1 = 1.0890
    
    --in nm
    Rm_S1 = 0.62193776E-09

    b_S1 = 0.06
    a_S1 = {-25933.587, 14.66188573699344914, 0.525743927693154455E06, -0.122790966318838728E07, 0.175565797136193828E06, 0.173795490253058379E07, -0.119112720845007316E07, -0.245659148870101490E07, 0.303380094883701415E08, -0.100054913157079869E09, -0.296340813141656632E08, 0.997302450614721887E09, -0.272673123492070958E10, 0.323269132716538832E10, -0.147953587185832486E10}
	
    
    r = r * rs
    ret = 0
    if r < R_SR_S1 then
        ret = h*c0 * (A_SR_S1 + B_SR_S1 * (a0 / r)^N_S1)
    elseif r >= R_SR_S1 and r <= R_LR_S1 then
        xi = (r - Rm_S1) / (r + b_S1 * Rm_S1)
        for i=0,14 do
            --print(i)
            ret = ret + a_S1[i+1] * xi^i
        end
        ret = ret * h * c0
    elseif r > R_LR_S1 then
        Vexch = h * c0 * Aex_S1 * (r / a0)^gamma_S1 * math.exp(- beta_S1 * r / a0)
        ret = -C6_S1 / (r^6) - C8_S1 / (r^8) - C10_S1 / (r^10) + Vexch
    end
    return ret
end
