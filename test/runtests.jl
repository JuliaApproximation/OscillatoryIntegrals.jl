using ApproxFun, OscillatoryIntegrals, Base.Test

f=Fun(exp)
x=Fun(identity)


for ω in (0.,0.001,0.1,1.,10.5,100.,1000.5)
   @test norm(diff(fourierintegral(f,ω)*exp(im*ω*x)) -f*exp(im*ω*x)) < 1000eps()
   @test_approx_eq sum(f*exp(im*ω*x)) fourier(f,ω)
end
