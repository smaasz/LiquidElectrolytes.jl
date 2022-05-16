module bvtest


using PyPlot


function main(;
              R0=1.0,
              Δg=0.0,
              n=1,
              A0=1.0,
              R=1.0,
              T=1.0,
              F=1.0,
              V=10.0,
              β=0.5
              )

    A(U)=1.0/(R*T)*(Δg + A0 + n*F*(U))

    r(U)=R0*( exp(-β*A(U)/2) - exp((1-β)*A(U)/2) )

    
    U=-V:V/100:V
    plot(U,r.(U))
    grid(true)
             
end


end
