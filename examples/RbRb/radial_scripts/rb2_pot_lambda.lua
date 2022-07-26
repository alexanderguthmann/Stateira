function radial_func(r)
    rs = 1.0E-10
    a0 = 5.29177210903E-11
    Eh = 4.3597447222071E-18

    r = r * rs / a0
    ret = -1.002320649155/(r^3)+9.18291680*math.exp(-0.7196 * r)
    ret = 4.6E-06 * Eh *ret
    
    return ret
end

--res = radial_func(5)
--print(res)
