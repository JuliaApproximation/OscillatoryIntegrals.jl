module OscillatoryIntegrals
    using Base,ApproxFun


export fourier,fourierintegral


## fourierintegral calculates a function u so that u(x) exp(i ω x)
#  is an indefinite integral \int f(x) exp(i ω x) dx
#  the normalization is not defined
# This implementation is based on [Keller 1999]

#TODO: choose rounding using domain(f)
fourierintegral(f::Fun{Chebyshev},ω)=ω==0?integrate(f):
        [BasisFunctional(ceil(Integer,ω));Derivative()+im*ω]\[0.,f]

## fourier returns

function fourier(f,ω)
    u=fourierintegral(f,ω)
    last(u)*exp(im*ω)-first(u)*exp(-im*ω)
end

end #module
