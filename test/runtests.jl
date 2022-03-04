using Test
using LiquidElectrolytes
using Unitful,PhysicalConstants
import LiquidElectrolytes: @siunits, @phconstants

@phconstants AvogadroConstant
@siunits km
@testset "units" begin
    @test    AvogadroConstant==6.02214076e23
    @test    km==1000
end

