using DataFrames

function molar_mass(MI)
	sql = "	SELECT 	molar_mass
			FROM	isotopologues
			WHERE 	(molecule_id, local_id)
			IN 		(VALUES " * join(["(" * join("?"^length(t) , ',') * ")" for t in MI], ',') * ")"			
	result = query_local_db(sql, [i for t in MI for i in t])
	
	out = []
	for row in result
		push!(out, row.molar_mass)
	end
	
	return out
end

## Default data for molecules and isotopologues for new databases
const molecules = DataFrame(id=Int[], formula=String[], name=String[])
push!(molecules,
    (1, "H2O", "Water"),
    (2, "CO2", "Carbon Dioxide"),
    (3,	"O3", "Ozone"),
    (4,	"N2O", "Nitrogen oxide"),
    (5,	"CO", "Carbon Monoxide"),
    (6,	"CH4", "Methane"),
    (7, "O2", "Oxygen"),
    (8,	"NO", "Nitric Oxide"),
    (9,	"SO2", "Sulfur Dioxide"),
    (10, "NO2", "Nitrogen Dioxide"),
    (11, "NH3", "Ammonia"),
    (12, "HNO3", "Nitric Acid"),
    (13, "OH", "Hydroxyl"),
    (14, "HF", "Hydrogen Fluoride"),
    (15, "HCl", "Hydrogen Chloride"),
    (16, "HBr", "Hydrogen Bromide"),
    (17, "HI", "Hydrogen Iodide"),
    (18, "ClO", "Chlorine Monoxide"),
    (19, "OCS", "Carbonyl Sulfide"),
    (20, "H2CO", "Formaldehyde"),
    (21, "HOCl", "Hypochlorous Acid"),
    (22, "N2", "Nitrogen"),
    (23, "HCN", "Hydrogen Cyanide"),
    (24, "CH3Cl", "Methyl Chloride"),
    (25, "H2O2", "Hydrogen Peroxide"),
    (26, "C2H2", "Acetylene"),
    (27, "C2H6", "Ethane"),
    (28, "PH3", "Phosphine"),
    (29, "COF2", "Carbonyl Fluoride"),
    (30, "SF6", "Sulfur Hexafluoride"),
    (31, "H2S", "Hydrogen Sulfide"),
    (32, "HCOOH", "Formic Acid"),
    (33, "HO2", "Hydroperoxyl"),
    (34, "O", "Oxygen Atom"),
    (35, "ClONO2", "Chlorine Nitrate"),
    (36, "NO+", "Nitric Oxide Cation"),
    (37, "HOBr", "Hypobromous Acid"),
    (38, "C2H4", "Ethylene"),
    (39, "CH3OH", "Methanol"),
    (40, "CH3Br", "Methyl Bromide"),
    (41, "CH3CN", "Acetonitrile"),
    (42, "CF4", "PFC-14"),
    (43, "C4H2", "Diacetylene"),
    (44, "HC3N", "Cyanoacetylene"),
    (45, "H2", "Hydrogen"),
    (46, "CS", "Carbon Monosulfide"),
    (47, "SO3", "Sulfur trioxide"),
    (48, "C2N2", "Cyanogen"),
    (49, "COCl2", "Phosgene"),
    (53, "CS2", "Carbon disulfide")
)

const isotopologues = DataFrame(
    molecule_id=Int[],
    global_id=Int[],
    local_id=Int[],
    formula=String[],
    afgl_code=Int[],
    abundance=Float64[],
    molar_mass=Float64[],
    q_t0=Float64[],
    g_i=Int[]
    )

push!(isotopologues,
    # 1: H2O
	(1, 1, 1, "H216O", 161, 0.997317, 18.010565, 174.58, 1),
	(1, 2, 2, "H218O", 181, 0.002000, 20.014811, 176.05, 1),
	(1, 3, 3, "H217O", 171, 3.718840e-4, 19.01478, 1052.14, 6),
	(1, 4, 4, "HD16O", 162, 3.106930e-4, 19.01674, 864.74, 6),
	(1, 5, 5, "HD18O", 182, 6.230030e-7, 21.020985, 875.57, 6),
	(1, 6, 6, "HD17O", 172, 1.158530e-7, 20.020956, 5226.79, 36),
	(1, 129, 7, "D216O", 262, 2.419700e-8, 20.022915, 1027.80, 1),
	
    # 2: CO2
	(2, 7, 1, "12C16O2", 626, 0.984204, 43.98983, 286.09, 1),
	(2, 8, 2, "13C16O2", 636, 0.011057, 44.993185, 576.64, 2),
	(2, 9, 3, "16O12C18O", 628, 0.003947, 45.994076, 607.81, 1),
	(2, 10, 4, "16O12C17O", 627, 7.339890e-4, 44.994045, 3542.61, 6),
	(2, 11, 5, "16O13C18O", 638, 4.434460e-5, 46.997431, 1225.46, 2),
	(2, 12, 6, "16O13C17O", 637, 8.246230e-6, 45.9974, 7141.32, 12),
	(2, 13, 7, "12C18O2", 828, 3.957340e-6, 47.998322, 323.42, 1),
	(2, 14, 8, "17O12C18O", 827, 1.471800e-6, 46.998291, 3766.58, 6),
	(2, 121, 9, "12C17O2", 727, 1.368470e-7, 45.998262, 10971.57, 1),
	(2, 15, 0, "13C18O2", 838, 4.446000e-8, 49.001675, 652.24, 2),
	(2, 120, 11, "18O13C17O", 837, 1.653540e-8, 48.001646, 7595.04, 12),
	(2, 122, 12, "13C17O2", 737, 1.537500e-9, 47.0016182378, 22120.47, 2),
	
    # 3: O3
	(3, 16, 1, "16O3", 666, 0.992901, 47.984745, 3483.71, 1),
	(3, 17, 2, "16O16O18O", 668, 0.003982, 49.988991, 7465.68, 1),
	(3, 18, 3, "16O18O16O", 686, 0.001991, 49.988991, 3647.08, 1),
	(3, 19, 4, "16O16O17O", 667, 7.404750e-4, 48.98896, 43330.85, 6),
	(3, 20, 5, "16O17O16O", 676, 3.702370e-4, 48.98896, 21404.96, 6),
	
    # 4: N2O
	(4, 21, 1, "14N216O", 446, 0.990333, 44.001062, 4984.90, 9),
	(4, 22, 2, "14N15N16O", 456, 0.003641, 44.998096, 3362.01, 6),
	(4, 23, 3, "15N14N16O", 546, 0.003641, 44.998096, 3458.58, 6),
	(4, 24, 4, "14N218O", 448, 0.001986, 46.005308, 5314.74, 9),
	(4, 25, 5, "14N217O", 447, 3.692800e-4, 45.005278, 30971.79, 54),
	
    # 5: CO
	(5, 26, 1, "12C16O", 26, 0.986544, 27.994915, 107.42, 1),
	(5, 27, 2, "13C16O", 36, 0.011084, 28.99827, 224.69, 2),
	(5, 28, 3, "12C18O", 28, 0.001978, 29.999161, 112.77, 1),
	(5, 29, 4, "12C17O", 27, 3.678670e-4, 28.99913, 661.17, 6),
	(5, 30, 5, "13C18O", 38, 2.222500e-5, 31.002516, 236.44, 2),
	(5, 31, 6, "13C17O", 37, 4.132920e-6, 30.002485, 1384.66, 12),
	
    # 6: CH4
	(6, 32, 1, "12CH4", 211, 0.988274, 16.0313, 590.48, 1),
	(6, 33, 2, "13CH4", 311, 0.011103, 17.034655, 1180.82, 2),
	(6, 34, 3, "12CH3D", 212, 6.157510e-4, 17.037475, 4794.73, 3),
	(6, 35, 4, "13CH3D", 312, 6.917850e-6, 18.04083, 9599.16, 6),
	
    # 7: O2
	(7, 36, 1, "16O2", 66, 0.995262, 31.98983, 215.73, 1),
	(7, 37, 2, "16O18O", 68, 0.003991, 33.994076, 455.23, 1),
	(7, 38, 3, "16O17O", 67, 7.422350e-4, 32.994045, 2658.12, 6),
	
    # 8: NO
	(8, 39, 1, "14N16O", 46, 0.993974, 29.997989, 1142.13, 3),
	(8, 40, 2, "15N16O", 56, 0.003654, 30.995023, 789.26, 2),
	(8, 41, 3, "14N18O", 48, 0.001993, 32.002234, 1204.44, 3),
	
    # 9: SO2
	(9, 42, 1, "32S16O2", 626, 0.945678, 63.961901, 6340.30, 1),
	(9, 43, 2, "34S16O2", 646, 0.041950, 65.957695, 6368.98, 1),
	
    # 10: NO2
	(10, 44, 1, "14N16O2", 646, 0.991616, 45.992904, 13577.48, 3),
	(10, 130, 2, "15N16O2", 656, 0.003646, 46.989938, 9324.70, 2),
	
    # 11: NH3
	(11, 45, 1, "14NH3", 4111, 0.995872, 17.026549, 1725.22, 3),
	(11, 46, 2, "15NH3", 5111, 0.003661, 18.023583, 1153.30, 2),
	
    # 12: HNO3
	(12, 47, 1, "H14N16O3", 146, 0.989110, 62.995644, 2.14e5, 6),
	(12, 117, 2, "H15N16O3", 156, 0.003636, 63.99268, 1.43e5, 4),
	
    # 13: OH
	(13, 48, 1, "16OH", 61, 0.997473, 17.00274, 80.35, 2),
	(13, 49, 2, "18OH", 81, 0.002000, 19.006986, 80.88, 2),
	(13, 50, 3, "16OD", 62, 1.553710e-4, 18.008915, 209.32, 3),
	
    # 14: HF
	(14, 51, 1, "H19F", 19, 0.999844, 20.006229, 41.47, 4),
	(14, 110, 2, "D19F", 29, 1.557410e-4, 21.012404, 115.91, 6),
	
    # 15: HCl
	(15, 52, 1, "H35Cl", 15, 0.757587, 35.976678, 160.65, 8),
	(15, 53, 2, "H37Cl", 17, 0.242257, 37.973729, 160.89, 8),
	(15, 107, 3, "D35Cl", 25, 1.180050e-4, 36.982853, 462.78, 12),
	(15, 108, 4, "D37Cl", 27, 3.773500e-5, 38.979904, 464.13, 12),
    
    # 16: HBr
	(16, 54, 1, "H79Br", 19, 0.506781, 79.92616, 200.17, 8),
	(16, 55, 2, "H81Br", 11, 0.493063, 81.924115, 200.23, 8),
	(16, 111, 3, "D79Br", 29, 7.893840e-5, 80.932336, 586.40, 12),
	(16, 112, 4, "D81Br", 21, 7.680160e-5, 82.930289, 586.76, 12),
	
    # 17: HI
	(17, 56, 1, "H127I", 17, 0.999844, 127.912297, 388.99, 12),
	(17, 113, 2, "D127I", 27, 1.557410e-4, 128.918472, 1147.06, 18),
	
    # 18: ClO
	(18, 57, 1, "35Cl16O", 56, 0.755908, 50.963768, 3274.61, 4),
	(18, 58, 2, "37Cl16O", 76, 0.241720, 52.960819, 3332.29, 4),
	
    # 19: OCS
	(19, 59, 1, "16O12C32S", 622, 0.937395, 59.966986, 1221.01, 1),
	(19, 60, 2, "16O12C34S", 624, 0.041583, 61.96278, 1253.48, 1),
	(19, 61, 3, "16O13C32S", 632, 0.010531, 60.970341, 2484.15, 2),
	(19, 62, 4, "16O12C33S", 623, 0.007399, 60.966371, 4950.11, 4),
	(19, 63, 5, "18O12C32S", 822, 0.001880, 61.971231, 1313.78, 1),
	(19, 135, 6, "16O13C34S", 634, 4.675080e-4, 62.966136, 2546.53, 2),
	
    # 20: H2CO
	(20, 64, 1, "H212C16O", 126, 0.986237, 30.010565, 2844.53, 1),
	(20, 65, 2, "H213C16O", 136, 0.011080, 31.01392, 5837.69, 2),
	(20, 66, 3, "H212C18O", 128, 0.001978, 32.014811, 2986.44, 1),
	
    # 21: HOCl
	(21, 67, 1, "H16O35Cl", 165, 0.755790, 51.971593, 19274.79, 8),
	(21, 68, 2, "H16O37Cl", 167, 0.241683, 53.968644, 19616.20, 8),
	
    # 22: N2
	(22, 69, 1, "14N2", 44, 0.992687, 28.006148, 467.10, 1),
	(22, 118, 2, "14N15N", 45, 0.007478, 29.003182, 644.10, 6),
	
    # 23: HCN
	(23, 70, 1, "H12C14N", 124, 0.985114, 27.010899, 892.20, 6),
	(23, 71, 2, "H13C14N", 134, 0.011068, 28.014254, 1830.97, 12),
	(23, 72, 3, "H12C15N", 125, 0.003622, 28.007933, 615.28, 4),
	
    # 24: CH3Cl
	(24, 73, 1, "12CH335Cl", 215, 0.748937, 49.992328, 57916.12, 4),
	(24, 74, 2, "12CH337Cl", 217, 0.239491, 51.989379, 58833.90, 4),
	
    # 25: H2O2
	(25, 75, 1, "H216O2", 1661, 0.994952, 34.00548, 9847.99, 1),
	
    # 26: C2H2
	(26, 76, 1, "12C2H2", 1221, 0.977599, 26.01565, 412.45, 1),
	(26, 77, 2, "H12C13CH", 1231, 0.021966, 27.019005, 1656.18, 8),
	(26, 105, 3, "H12C12CD", 1222, 3.045500e-4, 27.021825, 1581.84, 6),
    
    # 27: C2H6
	(27, 78, 1, "12C2H6", 1221, 0.976990, 30.04695, 70882.52, 1),
	(27, 106, 2, "12CH313CH3", 1231, 0.021953, 31.050305, 36191.80, 2),
	
    # 28: PH3
	(28, 79, 1, "31PH3", 1111, 0.999533, 33.997238, 3249.44, 2),
	
    # 29: COF2
	(29, 80, 1, "12C16O19F2", 269, 0.986544, 65.991722, 70028.43, 1),
	(29, 119, 2, "13C16O19F2", 369, 0.011083, 66.995083, 1.40e5, 2),
	
    # 30: SF6
	(30, 126, 1, "32S19F6", 29, 0.950180, 145.962492, 1.62e6, 1),
	
    # 31: H2S
	(31, 81, 1, "H232S", 121, 0.949884, 33.987721, 505.79, 1),
	(31, 82, 2, "H234S", 141, 0.042137, 35.983515, 504.35, 1),
	(31, 83, 3, "H233S", 131, 0.007498, 34.987105, 2014.94, 4),
	
    # 32: HCOOH
	(32, 84, 1, "H12C16O16OH", 126, 0.983898, 46.00548, 39132.76, 4),
	
    # 33: HO2
	(33, 85, 1, "H16O2", 166, 0.995107, 32.997655, 4300.39, 2),
	
    # 34: O
	(34, 86, 1, "16O", 6, 0.997628, 15.994915, 6.72, 1),
	
    # 35: ClONO2
	(35, 127, 1, "35Cl16O14N16O2", 5646, 0.749570, 96.956672, 4.79e6, 12),
	(35, 128, 2, "37Cl16O14N16O2", 7646, 0.239694, 98.953723, 4.91e6, 12),
	
    # 36: NO+
    (36, 87, 1, "14N16O+", 46, 0.993974, 29.997989, 311.69, 3),
	
    # 37: HOBr
	(37, 88, 1, "H16O79Br", 169, 0.505579, 95.921076, 28339.38, 8),
	(37, 89, 2, "H16O81Br", 161, 0.491894, 97.919027, 28237.98, 8),
	
    # 38: C2H4
	(38, 90, 1, "12C2H4", 221, 0.977294, 28.0313, 11041.54, 1),
	(38, 91, 2, "12CH213CH2", 231, 0.021959, 29.034655, 45196.89, 2),
	
    # 39: CH3OH
	(39, 92, 1, "12CH316OH", 2161, 0.985930, 32.026215, 70569.92, 2),
	
    # 40: CH3Br
	(40, 93, 1, "12CH379Br", 219, 0.500995, 93.941811, 83051.98, 4),
	(40, 94, 2, "12CH381Br", 211, 0.487433, 95.939764, 83395.21, 4),
	
    # 41: CH3CN
	(41, 95, 1, "12CH312C14N", 2124, 0.973866, 41.026549, 88672.19, 3),
	
    # 42: CF4
	(42, 96, 1, "12C19F4", 29, 0.988890, 87.993616, 1.21e5, 1),
	
    # 43: C4H2
	(43, 116, 1, "12C4H2", 2211, 0.955998, 50.01565, 9818.97, 1),
	
    # 44: HC3N
	(44, 109, 1, "H12C314N", 1224, 0.963346, 51.010899, 24786.84, 6),
	
    # 45: H2
	(45, 103, 1, "H2", 11, 0.999688, 2.01565, 7.67, 1),
	(45, 115, 2, "HD", 12, 3.114320e-4, 3.021825, 29.87, 6),
	
    # 46: CS
	(46, 97, 1, "12C32S", 22, 0.939624, 43.971036, 253.62, 1),
	(46, 98, 2, "12C34S", 24, 0.041682, 45.966787, 257.77, 1),
	(46, 99, 3, "13C32S", 32, 0.010556, 44.974368, 537.50, 2),
	(46, 100, 4, "12C33S", 23, 0.007417, 44.970399, 1022.97, 4),
	
    # 47: SO3
	(47, 114, 1, "32S16O3", 26, 0.943400, 79.95682, 7783.30, 1),
	
    # 48: C2N2
	(48, 123, 1, "12C214N2", 4224, 0.970752, 52.006148, 15582.44, 1),
	
    # 49: COCl2
	(49, 124, 1, "12C16O35Cl2", 2655, 0.566392, 97.9326199796, 1.48e6, 1),
	(49, 125, 2, "12C16O35Cl37Cl", 2657, 0.362235, 99.9296698896, 3.04e6, 16),
	
    # 53: CS2
	(53, 131, 1, "12C32S2", 222, 0.892811, 75.94414, 1352.60, 1),
	(53, 132, 2, "32S12C34S", 224, 0.079260, 77.93994, 2798.00, 1),
	(53, 133, 3, "32S12C33S", 223, 0.014094, 76.943256, 1107.00, 4),
	(53, 134, 4, "13C32S2", 232, 0.010310, 76.947495, 2739.70, 2),
)
