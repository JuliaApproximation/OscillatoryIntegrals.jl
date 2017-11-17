__precompile__()
module OscillatoryIntegrals
    using Base, ApproxFun

    import ApproxFun: domain, evaluate, spacescompatible,
                        SpaceOperator, ConstantSpace

    include("Bessel.jl")

# TODO: Fourier weight

export fourier,fourierintegral,fouriercauchy


## fourierintegral calculates a function u so that u(x) exp(i ω x)
#  is an indefinite integral \int f(x) exp(i ω x) dx
#  the normalization is not defined
# This implementation is based on [Keller 1999]

#TODO: choose rounding using domain(f)
function fourierintegral(f::Fun{<:Chebyshev}, ω)
    if ω == 0
        integrate(f)
    else
        \([BasisFunctional(ceil(Integer,ω),space(f));
         Derivative(space(f))+im*ω*I], [0.,f]; tolerance=10maximum(f.coefficients)*eps())
    end
end

## fourier returns the fourier transform
# ∫ f(x) exp(i*w*x) dx

function fourier(sp::Space, f, ω)
    d = domain(sp)
    u = fourierintegral(Fun(sp,f),ω)
    last(u)*exp(im*ω*last(d))-first(u)*exp(im*ω*first(d))
end

# Fourier transform of Legendre polynomials is known
# ∫₋₁⁺¹ Pⱼ(x) e⁻ⁱʷˣ dx = 2(-i)ʲ jⱼ(ω)

function fourier(f::Fun{Jacobi{Segment{T},T},T}) where T
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
    Fun(Besselj(),cfs)
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


function fourier(S::Laurent{DD,RR}, f, ω) where {DD<:PeriodicLine,RR}
    #@assert domain(f)==PeriodicLine()
    if ω>0
        4π*exp(-ω)*webersum(f[2:2:end],ω)
    else # ω<0
        4π*exp(ω)*webersum(f[3:2:end],-ω)
    end
end

function fourier(S::Space{<:Ray}, f, ω)
    D=Derivative(S)
    u=(D+im*ω)\Fun(S,f)
    -first(u)
end

function fourier(S::Space{<:Line},f,ω)
    #@assert domain(f)==Line()
    fourier(Fun(x->evaluate(f,S,x),Ray()),ω) -
        fourier(Fun(x->evaluate(f,S,x),Ray(0.,π)),ω)
end



## fouriercauchy
# calculates cauchy transform of f(x)*exp(im*ω*x)

fouriercauchy(s::Bool,f::Fun,ω,z) = fouriercauchy(s,space(f),coefficients(f),ω,z)

function fouriercauchy(s::Bool,S::Space{<:Line},f,ω,z)
    @assert isreal(z)
    x=Fun(identity,S)
    M=Multiplication(x-z,S)
    g=Fun(rangespace(M), f.-evaluate(f,S,z))
    u=M\g
    ret=fourier(u,ω)/(2π*im)
    s ? (ret+evaluate(f,S,z)*exp(im*ω*z)) : ret
end


end #module
