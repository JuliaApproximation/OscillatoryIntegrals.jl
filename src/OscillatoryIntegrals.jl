module OscillatoryIntegrals
    using Base,ApproxFun


# TODO: Fourier weight

export fourier,fourierintegral,fouriercauchy


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
    p0=zero(ω)
    p1=one(ω)
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


function fourier{S,R<:Ray,T}(f::Fun{MappedSpace{S,R,T}},ω)
    if ω==0
        sum(f)
    else
        D=Derivative(space(f))
        u=(D+im*ω)\f
        -first(u)
    end
end

function fourier{S,L<:Line,T}(f::Fun{MappedSpace{S,L,T}},ω)
    @assert domain(f)==Line()
    if ω==0
        sum(f)
    else
        fourier(Fun(x->f[x],Ray()),ω)-fourier(Fun(x->f[x],Ray(0.,π)),ω)
    end
end



## fouriercauchy
# calculates cauchy transform of f[x]*exp(im*ω*x)


function fouriercauchy{S,L<:Line,T}(s::Bool,f::Fun{MappedSpace{S,L,T}},ω,z)
    @assert isreal(z)
    x=Fun(identity,space(f))
    M=Multiplication(x-z,space(f))
    g=Fun(f-f[z],rangespace(M))
    u=M\g
    ret=fourier(u,ω)/(2π*im)
    s?(ret+f[z]*exp(im*ω*z)):ret
end


end #module
