function radial_func(r)
    De = 8516.709
    re = 2.672993
    c_arr ={6.71527E06, 1.12588E08, 2.78604E09}
    p = 5
    q = 3
    rref = 3.85
    b_arr = {-2.89828701, -1.309265, -2.018507, -1.38162, -1.21933, 0.3463, 0.1061, -0.1886, 1.671, 10.943, 3.944, -27.23, -11.49, 56.7, 38.9, -37.4, -33.0}

    pad = 6
    qad = 6

    u_infty = 0
    u_arr = {-0.03228048156080912, 0.0016639423484953154, -0.0648937515913173}


    S = 3.818471145493912
    
    --Regular MLR Potential:
        
    yp = (r^p - rref^p) / (r^p + rref^p)
    yq = (r^q - rref^q) / (r^q + rref^q)
    yp_eq = (r^p - re^p) / (r^p + re^p)
    
    u_lr_r = 0
    u_lr_re = 0
    for i=1,3 do
        --print(i)
        u_lr_r = u_lr_r + c_arr[i] / (r^(6+2*(i-1)))
        u_lr_re = u_lr_re + c_arr[i] / (re^(6+2*(i-1)))
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

--res = radial_func(10)
--print(res)
