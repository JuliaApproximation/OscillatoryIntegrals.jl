using ApproxFun, OscillatoryIntegrals, Base.Test

f=Fun(exp)
x=Fun(identity)


for ω in (0.,0.001,0.1,1.,10.5,100.,1000.5)
   @test norm(diff(fourierintegral(f,ω)*exp(im*ω*x)) -f*exp(im*ω*x)) < 1000eps()
   @test_approx_eq sum(f*exp(im*ω*x)) fourier(f,ω)
end


f=Fun(sech,Laurent(PeriodicLine()))
@test_approx_eq_eps fourier(f,10.4) 5.051528742182663e-7 10e-13
@test_approx_eq_eps fourier(f,-10.4) 5.051528742182663e-7 10e-13

f=Fun(x->sech(x-0.1),Laurent(PeriodicLine()))
@test_approx_eq_eps fourier(f,10.4) (2.557186179286481e-7 + 4.356459741299553e-7*im) 10e-13
@test_approx_eq_eps fourier(f,-10.4) (2.557186179286481e-7 - 4.356459741299553e-7*im) 10e-13

