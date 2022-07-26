function radial_func(r)
    De = 333.758
    re = 4.17005
    c_arr ={6718500, 112629000, 2786830000}
    p = 5
    q = 3
    rref = 8.0
    rho = 0.54
    b_arr = {-0.51608, -0.09782, 0.1137, -0.0248}

    pad = 6
    qad = 6

    u_infty = 0
    u_arr = {-0.00981725985612236}

    S = 1.5389576619468355
    
    --Regular MLR Potential:
        
    yp = (r^p - rref^p) / (r^p + rref^p)
    yq = (r^q - rref^q) / (r^q + rref^q)
    yp_eq = (r^p - re^p) / (r^p + re^p)

    bds = 3.30
    cds = 0.423
    
    u_lr_r = 0
    u_lr_re = 0
    for i=1,3 do
        --print(i)
        dds_re = (1.0 - math.exp(-1.0*( bds*rho*re / (6+2*(i-1)) + cds * (rho*re)^2 / math.sqrt(6+2*(i-1)) )))^(5+2*(i-1))
        dds_r = (1.0 - math.exp(-1.0*( bds*rho*r / (6+2*(i-1)) + cds * (rho*r)^2 / math.sqrt(6+2*(i-1)) )))^(5+2*(i-1))
        
        u_lr_r = u_lr_r + dds_r*c_arr[i] / (r^(6+2*(i-1)))
        u_lr_re = u_lr_re + dds_re*c_arr[i] / (re^(6+2*(i-1)))
    end
    
    beta_infty = math.log(2*De / u_lr_re)
    
    beta = 0
    for i=#b_arr,1,-1 do
        --print(i)
        beta = beta*yq + b_arr[i]
    end
    beta = beta_infty * yp + (1 - yp) * beta
    
    vv = De*(1.0 - u_lr_r / u_lr_re * math.exp(-1.0 * beta * yp_eq))^2 - De
    
    --BOB Correction term:
    ypad = (r^pad - re^pad) / (r^pad + re^pad)
    yqad = (r^qad - re^qad) / (r^qad + re^qad)
    
    yqadsm = 1.0 - ypad
    
    sc = u_infty * ypad
    
    for i=1,#u_arr do
        --print(i)
        sc = sc + yqadsm * u_arr[i]
        yqadsm = yqadsm * yqad
    end
        
    vv = vv + 2 *sc
    Ecm = 1.9864455E-23

    if r < re then
	Eh = 4.3597447222071E-18
	a0 = 5.29177210903E-01
	v_shift = S * (r - re)^2
	v_shift = v_shift * 1.0E-06 * Eh /(a0^2)
        
    else
	v_shift = 0
    end
        
    v_shift = v_shift

    return vv * Ecm + v_shift
end

--res = radial_func(8)
--print(res)
