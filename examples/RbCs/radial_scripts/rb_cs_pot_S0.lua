function radial_func(r)
    a0 = 5.29177210903E-11
    rs = 1.0E-10
    Eh = 4.3597447222071E-18
    h = 6.62607015E-34
    c0 = 2.99792458E08
    R_SR_S0 = 0.3315E-09
    R_LR_S0 = 1.150E-09
    A_SR_S0 = -0.407634031E06
    B_SR_S0 = 1.52770630E11
    N_S0 = 7
    --coefficients for long range potential
    C6_S0 = 5693.7056 * Eh * a0^6
    C8_S0 = 796487.36 * Eh * a0^8
    C10_S0 = 95332817 * Eh * a0^10
    Aex_S0 = 0.37664685E05
    gamma_S0 = 5.427916
    beta_S0 = 1.0890
    --in nm
    Rm_S0 = 0.442708150E-09


    b_S0 = 0.09
    a_S0 = {-383636.509, -3.69980716645394794, 0.447519742785341805E7, -0.134065881674135253E7, -0.112246913875781145E8, -0.680373468487243954E7, 0.124395856928352383E8, -0.527808915105630062E8, 0.160604050855185674E9, 0.856669313055434823E9, -0.423220682973604128E10, -0.846286860630152822E10, 0.775110557475278497E11, 0.208102060193851382E11, -0.762262944271048737E12, 0.645280096247728157E12, 0.358089708848128967E13, -0.685156406423631516E13, -0.340359743040435295E13, 0.204117122912590576E14, -0.207876500106921722E14, 0.712777331768994293E13}
	
    
    r = r * rs
    ret = 0
    if r < R_SR_S0 then
        ret = h*c0 * (A_SR_S0 + B_SR_S0 * (a0 / r)^N_S0)
    elseif r >= R_SR_S0 and r <= R_LR_S0 then
        xi = (r - Rm_S0) / (r + b_S0 * Rm_S0)
        for i=0,21 do
            --print(i)
            ret = ret + a_S0[i+1] * xi^i
        end
        ret = ret * h * c0
    elseif r > R_LR_S0 then
        Vexch = h * c0 * Aex_S0 * (r / a0)^gamma_S0 * math.exp(- beta_S0 * r / a0)
        ret = -C6_S0 / (r^6) - C8_S0 / (r^8) - C10_S0 / (r^10) - Vexch
    end
    return ret
end