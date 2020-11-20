using Test

cd("..")  # So that "deps/penelope.so" is the relative location of the object file.
using Penelope

@testset "load_so" begin
    @test Penelope.greet() == nothing
    E = 20000.0
    loc = [1,2,3]
    dir = [0,0,-1]
    Penelope.initialize_track(E, loc, dir)
    @test E == Penelope.energy()
    @test loc == Penelope.location()
    @test dir == Penelope.direction()
end
