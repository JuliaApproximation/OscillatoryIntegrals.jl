using ApproxFun, OscillatoryIntegrals, LinearAlgebra, Test

@testset "bounded interval" begin
   f=Fun(exp)
   x=Fun(identity)

   for ω in (0.,0.001,0.1,1.,10.5,100.,1000.5)
      @test norm((fourierintegral(f,ω)*exp(im*ω*x))' -f*exp(im*ω*x)) < 1000eps()
      @test sum(f*exp(im*ω*x)) ≈ fourier(f,ω)
   end

   x = Fun(0..1)
   @test fourier(x, 2.0) ≈ sum(x*exp(2im*x))
end

@testset "Periodic line" begin
   f=Fun(sech,Laurent(PeriodicLine()))
   @test fourier(f,10.4) ≈ 5.051528742182663e-7 atol=10e-13
   @test fourier(f,-10.4) ≈ 5.051528742182663e-7 atol=10e-13

   f=Fun(x->sech(x-0.1),Laurent(PeriodicLine()))
   @test fourier(f,10.4) ≈ (2.557186179286481e-7 + 4.356459741299553e-7*im) atol=10e-13
   @test fourier(f,-10.4) ≈ (2.557186179286481e-7 - 4.356459741299553e-7*im) atol=10e-13
end

# @test fourier(f,2.) ≈ 0.27101495139940088877738507029


@testset "Wiener–Hopf kernel" begin
   γ = 2.0
   Γ = Ray{true}(0.0,false) ∪ Ray() # Line split in two
   K = Fun(x -> exp(-γ*abs(x)), Γ)

   @test fourier(component(K,1), -1.0) ≈ 1/(γ-im*1.0)
   @test fourier(component(K,2), -1.0) ≈ 1/(γ+im*1.0)
   @test fourier(K, -1.0) ≈ 2γ/(γ^2+1.0^2)
   @test fourier(K, 2.0) ≈ 2γ/(γ^2+2.0^2)
end



