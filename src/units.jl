import Unitful
import PhysicalConstants

SI(x)=Float64(Unitful.ustrip(Unitful.upreferred(1*x)));

macro  si_str(x)
    Unitful=getproperty(@__MODULE__,:Unitful)
    quote
        Float64($(Unitful).ustrip($(Unitful).upreferred(1*$(Unitful).@u_str($(x)))))
    end
end

function _siunits(xs...)
    Unitful=getproperty(@__MODULE__,:Unitful)
    code = Expr(:block)
    for x in xs
        push!(code.args, :(const $x = Float64($(Unitful).ustrip($(Unitful).upreferred(1*$(Unitful).$x)))))
    end
    code
end

macro siunits(xs...)
    esc(_siunits(xs...))
end


function _phconstants(xs...)
    Unitful=getproperty(@__MODULE__,:Unitful)
    PhysicalConstants=getproperty(@__MODULE__,:PhysicalConstants)
    code = Expr(:block)
    for x in xs
        push!(code.args, :(const $x = Float64($(Unitful).ustrip($(Unitful).upreferred($(PhysicalConstants).CODATA2018.$x)))))
    end
    code
end

macro phconstants(xs...)
    esc(_phconstants(xs...))
end

