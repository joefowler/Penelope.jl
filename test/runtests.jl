using Test

cd("..")  # So that "deps/penelope.so" is the relative location of the object file.
using Penelope

@testset "load_so" begin
    @test Penelope.greet() == nothing
    want = 3
    unsafe_store!(Penelope.track.ptr_kpar, want)
    @test want == unsafe_load(Penelope.track.ptr_kpar)
end
