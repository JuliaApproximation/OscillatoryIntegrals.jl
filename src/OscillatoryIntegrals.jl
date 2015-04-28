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

## fourier returns the fourier transform
# \int f(x) exp(i*w*x) dx

function fourier(f,ω)
    u=fourierintegral(f,ω)
    last(u)*exp(im*ω)-first(u)*exp(-im*ω)
end



## Webers algorithm for Fourier transforms

function webersum(cfs,ω)
    p0=0
    p1=1
    ret=p1*first(cfs)
    for k=2:length(cfs)
        p1,p0=2*(1-k+ω)/(k-1)*p1-p0,p1
        ret+=p1*cfs[k]
    end
    ret
end

function fourier(f::Fun{Laurent},ω)
    @assert domain(f)==PeriodicLine()

    if ω==0
        sum(f)
    elseif ω>0
        4π*exp(-ω)*webersum(f.coefficients[2:2:end],ω)
    else # ω<0
        4π*exp(ω)*webersum(f.coefficients[3:2:end],-ω)
    end
end

end #module
