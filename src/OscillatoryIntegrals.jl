__precompile__()
module OscillatoryIntegrals
    using Base, ApproxFun, Plots

    import ApproxFun: UnivariateSpace, RealUnivariateSpace, domain, evaluate

    include("Bessel.jl")

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
# ∫ f(x) exp(i*w*x) dx

function fourier(f,ω)
    u=fourierintegral(f,ω)
    last(u)*exp(im*ω)-first(u)*exp(-im*ω)
end

# Fourier transform of Legendre polynomials is known
# ∫₋₁⁺¹ Pⱼ(x) e⁻ⁱʷˣ dx = 2(-i)ʲ jⱼ(ω)

function fourier{T}(f::Fun{Jacobi{T,Interval{T}},T})
    a,b = f.space.a,f.space.b
    @assert a == b == 0
    @assert domain(f) == Interval{T}()
    cfs = complex(coefficients(f))
    s = 2
    for i=1:2:length(cfs)
        cfs[i] *= s
        s *= -1
    end
    s = -2
    for i=2:2:length(cfs)
        cfs[i] *= im*s
        s *= -1
    end
    Fun(cfs,Besselj())
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

fourier(f::Fun,ω) = ω == 0 ? sum(f) : fourier(space(f),coefficients(f),ω)


function fourier{DD<:PeriodicLine}(S::Laurent{DD},f,ω)
    #@assert domain(f)==PeriodicLine()
    if ω>0
        4π*exp(-ω)*webersum(f[2:2:end],ω)
    else # ω<0
        4π*exp(ω)*webersum(f[3:2:end],-ω)
    end
end

function fourier{T,R<:Ray}(S::UnivariateSpace{T,R},f,ω)
    D=Derivative(S)
    u=(D+im*ω)\Fun(f,S)
    -first(u)
end

function fourier{T,L<:Line}(S::UnivariateSpace{T,L},f,ω)
    #@assert domain(f)==Line()
    fourier(Fun(x->evaluate(f,S,x),Ray()),ω)-fourier(Fun(x->evaluate(f,S,x),Ray(0.,π)),ω)
end



## fouriercauchy
# calculates cauchy transform of f[x]*exp(im*ω*x)

fouriercauchy(s::Bool,f::Fun,ω,z) = fouriercauchy(s,space(f),coefficients(f),ω,z)

function fouriercauchy{T,L<:Line}(s::Bool,S::UnivariateSpace{T,L},f,ω,z)
    @assert isreal(z)
    x=Fun(identity,S)
    M=Multiplication(x-z,S)
    g=Fun(f-evaluate(f,S,z),rangespace(M))
    u=M\g
    ret=fourier(u,ω)/(2π*im)
    s?(ret+evaluate(f,S,z)*exp(im*ω*z)):ret
end


end #module
