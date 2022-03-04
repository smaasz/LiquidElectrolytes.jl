import Unitful
import PhysicalConstants

macro siunits(xs...)
    code = Expr(:block)
    for x in xs
        push!(code.args, :(const $x = Float64(Unitful.ustrip(Unitful.upreferred(1Unitful.$x)))))
    end
    esc(quote eval($code) end)
end

macro phconstants(xs...)
    code = Expr(:block)
    for x in xs
        push!(code.args, :(const $x = Float64(Unitful.ustrip(Unitful.upreferred(PhysicalConstants.CODATA2018.$x)))))
    end
    esc(quote eval($code) end)
end

