function radial_func(r)
    rs = 1.0E-10
    a0 = 5.29177210903E-11
    Eh = 4.3597447222071E-18
    alpha = 7.2973525693E-03
    A_short = -50.974
    A_long = -0.0525
    beta_short = 0.8
    beta_long = 0.28

    r = r * rs
    ret = Eh * alpha^2 * ( A_short * math.exp(-beta_short * r / a0) + A_long * math.exp(-beta_long * r / a0) + 1.0/((r/a0)^3) )
    
    return ret
end

