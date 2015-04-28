# OscillatoryIntegrals.jl
Calculate oscillatory integrals using Julia


`
using ApproxFun,OscillatoryIntegrals
fourier(Fun(x->cos(x^2),[1,2]),10000.5)`

calculates the integral of `cos(x^2)*exp(im*10000.5*x)` over `[1,2]`

`fourier(Fun(sech,Laurent([-Inf,Inf])),10.5)`

calcules the integral of `sech(x)*exp(im*10.5*x)` over the real line