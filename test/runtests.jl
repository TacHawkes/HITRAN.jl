using HITRAN
using Test

@testset "HITRAN.jl" begin
    @test_nowarn test_db = open_database("test.sqlite")
    @testset "Total internal partition sums" begin
        for i in HITRAN.isotopologues[:, :global_id]
            q_t0 = HITRAN.isotopologues[HITRAN.isotopologues.global_id .== i, :q_t0][1]            
            tolerance = 10.0^(ceil(Int, log10(q_t0))-3)
            @test isapprox(tips(i, HITRAN.c_T_ref), q_t0, atol=tolerance)
        end
    end    
end

