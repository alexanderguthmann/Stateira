function radial_func(r)
    rs = 1.0E-10
    a0 = 5.29177210903E-11
    Eh = 4.3597447222071E-18
    alpha = 7.2973525693E-03


    r = r * rs
    ret = Eh * alpha^2 *  1.0/((r/a0)^3) 
    
    return ret 
end

