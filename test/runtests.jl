using HITRAN
using Test

@testset "HITRAN.jl" begin
    @testset "Total internal partition sums" begin for i in HITRAN.isotopologues[:,
                                                                                 :global_id]
        q_t0 = HITRAN.isotopologues[HITRAN.isotopologues.global_id .== i, :q_t0][1]
        tolerance = 10.0^(ceil(Int, log10(q_t0)) - 3)
        @test isapprox(tips(i, HITRAN.c_T_ref), q_t0, atol = tolerance)
    end end
    @testset "Database and absorption calculation" begin
        open_database(joinpath(@__DIR__, "HITRAN.sqlite"))
        @test_nowarn fetch!("StdAtm", iso_id(7, 1), 13000, 13150, [:standard])
        @test isequal(iso_id(["H2O"])[1], 1)
        @test isequal(iso_id(7, 1)[1], 36)
        @test_nowarn α(["StdAtm"])
        @test_nowarn α(["StdAtm"], :voigt)
        @test_nowarn α(["StdAtm"], :sdvoigt)
        @test_nowarn α(["StdAtm"], :lorentz)
        @test_nowarn α(["StdAtm"], :gauss)
        ν, absorption = α(["StdAtm"])
        @test_nowarn apply_instrument_function(ν, absorption, :gaussian, 1.0, 0.1)
        @test isapprox(absorption_spectrum([0.0], 1)[1], 0.0)
        @test isapprox(transmittance_spectrum([0.0], 1)[1], 1.0)
        @test isapprox(optical_depth([0.0], 1)[1], 0.0)
    end

    @testset "Instrument functions" begin
        x = -10:0.001:10
        for (s, fn) in HITRAN.instrument_functions
            @test isapprox(sum(fn(x, 0.1) * 0.001), 1, atol = 1e-2)
        end
    end

    @testset "Environments" begin
        @test isapprox(HITRAN.p_s_h2o(273.15), 611.2911778902558)
        @test isapprox(HITRAN.kelvin_to_celsius(273.15), 0)
        @test isapprox(HITRAN.celsius_to_kelvin(0), 273.15)
        @test_nowarn moist_air(100)
        @test_nowarn moist_air(100, HITRAN.c_p_ref, 240.0)
    end

    @testset "Utility" begin
        @test isapprox(wavelength_to_wavenumber(1e-6), 10000.0)
        @test isapprox(frequency_to_wavenumber(30e9), 1, atol = 1e-3)
    end
end
