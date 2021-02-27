using HITRAN
using Test

@testset "HITRAN.jl" begin    
    @testset "Total internal partition sums" begin
        for i in HITRAN.isotopologues[:, :global_id]
            q_t0 = HITRAN.isotopologues[HITRAN.isotopologues.global_id .== i, :q_t0][1]            
            tolerance = 10.0^(ceil(Int, log10(q_t0))-3)
            @test isapprox(tips(i, HITRAN.c_T_ref), q_t0, atol=tolerance)
        end
    end
    @testset "Database and absorption calculation" begin
        @test_nowarn fetch!("StdAtm", iso_id(7, 1), 13000, 13150, [:standard])
        @test isequal(iso_id(["H2O"])[1], 1)
        @test isequal(iso_id(7, 1)[1], 36)
        @test_nowarn α(["StdAtm"])
        @test_nowarn α(["StdAtm"], :voigt)
        @test_nowarn α(["StdAtm"], :sdvoigt)
        @test_nowarn α(["StdAtm"], :lorentz)
        @test_nowarn α(["StdAtm"], :gauss)
        @test isapprox(absorption_spectrum([0.], 1)[1], 0.0)
        @test isapprox(transmittance_spectrum([0.], 1)[1], 1.0)
        @test isapprox(optical_depth([0.], 1)[1], 0.0)
    end

    @testset "Water saturation pressure" begin
        @test isapprox(HITRAN.p_s_h2o(273.15), 611.2911778902558)
    end
end

