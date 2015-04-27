# OscillatoryIntegrals.jl
Calculate oscillatory integrals using Julia


```
using ApproxFun,OscillatoryIntegrals
fourier(Fun(x->cos(x^2),[1,2]),10000.5)```

calculates the integral of `cos(x^2)*exp(im*10000.5*x)`
