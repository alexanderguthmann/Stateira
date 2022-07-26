function radial_func(r)
    a0 = 5.29177210903E-11
    rs = 1.0E-10
    Eh = 4.3597447222071E-18
    h = 6.62607015E-34
    c0 = 2.99792458E08
    R_SR_S1 = 5.07E-10
    R_LR_S1 = 11E-10
    N_S1 = 4.5338950
    A_SR_S1 = -0.619088543E05
    B_SR_S1 = 0.956231677E08 * (1.0E-10)^N_S1
    --coefficients for long range potential
    C6 = 0.2270032E10 * (1.0E-10)^6
    C8 = 0.7782886E11 * (1.0E-10)^8
    C10 = 0.2868869E13 * (1.0E-10)^10
    C26 = 0.2819810E28 * (1.0E-10)^26
    gamma = 5.317689
    beta = 2.093816 * (1.0E-10)^-1
    Aex = 0.1317785E07 * (1.0E-10)^(-gamma)



    
    --in nm
    Rm_S1 = 6.0933451E-10

    b_S1 = -0.33
    a_S1 = {-24150.3352, -67.2503402304666542, 0.195494577140503543E6, -0.141544168453406223E6, -0.221166468149940465E6, 0.165443726445793004E6, -0.596412188910614259E6, 0.654481694231538040E6, 0.261413416681972012E7, -0.349701859112702878E7, -0.328185277155018630E7, 0.790208849885562522E7, -0.398783520249289213E7}
    
	
    
    r = r * rs
    ret = 0
    if r < R_SR_S1 then
        ret = h*c0 * (A_SR_S1 + B_SR_S1 * (1.0 / r)^N_S1)
    elseif r >= R_SR_S1 and r <= R_LR_S1 then
        xi = (r - Rm_S1) / (r + b_S1 * Rm_S1)
        for i=0,12 do
            --print(i)
            ret = ret + a_S1[i+1] * xi^i
        end
        ret = ret * h * c0
    elseif r > R_LR_S1 then
        Vexch = Aex * r^gamma * math.exp(-beta * r )
        ret =  h * c0 * (-C6 / (r^6) - C8 / (r^8) - C10 / (r^10) - C26 / (r^26) + Vexch)
    end
    return ret
end

--res = radial_func(0.5)
--print(res)