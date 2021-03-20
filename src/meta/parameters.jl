function merge_groups(g::Symbol...)        
    return unique(vcat(getindex.(Ref(parameter_groups), g)...))
end

# initialize parameter groups with all basic sets and then add combined ones
const parameter_groups = Dict(
    :standard => ["molec_id", "local_iso_id", "nu", "sw", "a", "elower", "gamma_air",
                "delta_air","gamma_self","n_air","n_self","gp","gpp"],

    :labels => ["statep", "statepp"],

    :voigt_air => ["gamma_air","delta_air","deltap_air","n_air"],
    :voigt_self => ["gamma_self","delta_self","deltap_self","n_self"],
    :voigt_h2 => ["gamma_H2","delta_H2","deltap_H2","n_H2"],
    :voigt_co2 => ["gamma_CO2","delta_CO2","n_CO2"],
    :voigt_he => ["gamma_He","delta_He","n_He"],
    :voigt_h2o => ["gamma_H2O","n_H2O"],
    :voigt_linemixing => ["y_air","y_self"],

    :sdvoigt_air => ["gamma_air","delta_air","deltap_air","n_air", "SD_air"],
    :sdvoigt_self => ["gamma_self","delta_self","deltap_self","n_self", "SD_self"],    
    :sdvoigt_linemixing => ["Y_SDV_air_296","Y_SDV_self_296"],

    :ht_self => ["gamma_HT_0_self_50","n_HT_self_50","gamma_HT_2_self_50",
                "delta_HT_0_self_50","deltap_HT_self_50","delta_HT_2_self_50",
                "gamma_HT_0_self_150","n_HT_self_150","gamma_HT_2_self_150",
                "delta_HT_0_self_150","deltap_HT_self_150","delta_HT_2_self_150",
                "gamma_HT_0_self_296","n_HT_self_296","gamma_HT_2_self_296",
                "delta_HT_0_self_296","deltap_HT_self_296","delta_HT_2_self_296",
                "gamma_HT_0_self_700","n_HT_self_700","gamma_HT_2_self_700",
                "delta_HT_0_self_700","deltap_HT_self_700","delta_HT_2_self_700",
                "nu_HT_self","kappa_HT_self","eta_HT_self","Y_HT_self_296"],                
    :ht_air => ["gamma_HT_0_air_296","n_HT_air_296","gamma_HT_2_air_296",
                "delta_HT_0_air_296","deltap_HT_air_296","delta_HT_2_air_296",
                "nu_HT_air","kappa_HT_air","eta_HT_air","Y_HT_air_296"]

)
parameter_groups[:voigt_all] = merge_groups(:voigt_air, :voigt_self, :voigt_h2, :voigt_co2, :voigt_he, :voigt_h2o, :voigt_linemixing)
parameter_groups[:sdvoigt_all] = merge_groups(:sdvoigt_air, :sdvoigt_self, :sdvoigt_linemixing)
parameter_groups[:ht_all] = merge_groups(:ht_self, :ht_air)

parameter_groups[:all] = merge_groups(:standard, :labels, :voigt_all, :sdvoigt_all, :ht_all)