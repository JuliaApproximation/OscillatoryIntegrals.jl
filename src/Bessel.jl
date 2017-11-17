export Bessel
export BesselJ, BesselY, BesselI, BesselK
export Besselj, Bessely, Besseli, Besselk

struct Bessel{KIND,MOD,GEOM} <: Space{Line{false,Float64},Float64} end

const BesselJ = Bessel{1,false,:cylindrical}
const BesselY = Bessel{2,false,:cylindrical}
const BesselI = Bessel{1,true,:cylindrical}
const BesselK = Bessel{2,true,:cylindrical}
const Besselj = Bessel{1,false,:spherical}
const Bessely = Bessel{2,false,:spherical}
const Besseli = Bessel{1,true,:spherical}
const Besselk = Bessel{2,true,:spherical}

domain(B::Bessel) = ℝ
spacescompatible{KIND,MOD,GEOM}(A::Bessel{KIND,MOD,GEOM},B::Bessel{KIND,MOD,GEOM}) = true

for (KIND,MOD,JYIK,jyik) in ((1,false,"J","j"),(2,false,"Y","y"),(1,true,"I","i"),(2,true,"K","k"))
    bjyik = parse(":bessel$jyik")
    @eval begin
        Base.show(io::IO,B::Bessel{$KIND,$MOD,:cylindrical}) = print(io,string($JYIK,"ᵢ(ℝ)"))
        Base.show(io::IO,B::Bessel{$KIND,$MOD,:spherical}) = print(io,string($jyik,"ᵢ(ℝ)"))
        evaluate(f::AbstractVector,S::Bessel{$KIND,$MOD,:cylindrical},ω::Number) = ω ≥ 0 ? dotu(f,$bjyik(0:length(f)-1,ω)) : dotu(f,(-1.).^(0:length(f)-1).*$bjyik(0:length(f)-1,-ω))
        evaluate(f::AbstractVector,S::Bessel{$KIND,$MOD,:spherical},ω::Number) = ω == 0 ? first(f) : ω > 0 ? dotu(f,$bjyik(0.5:length(f)-0.5,ω)*sqrt(π/2ω)) : dotu(f,(-1.).^(0:length(f)-1).*$bjyik(0.5:length(f)-0.5,-ω)*sqrt(π/(-2ω)))
        evaluate(f::AbstractVector,S::Bessel{$KIND,$MOD},ω::AbstractVector) = [evaluate(f,S,ωk) for ωk in ω]
    end
end

# Plots.plot{B<:Bessel,T<:Real}(f::Fun{B,Complex{T}};grid=true,kwds...) = plot!(plot(grid=grid),f;kwds...)
# Plots.plot!{B<:Bessel,T<:Real}(plt::Plots.Plot,f::Fun{B,Complex{T}};kwds...) = plot!(plt,linspace(-max(length(f),10),max(length(f),10),6length(f)+100),f;kwds...)
