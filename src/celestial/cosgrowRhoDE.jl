# Date created: 10/01/2020
# Julia Conversion: Gerry Gralton
# Original author: Aaron Robotham

function cosgrowRhoDE(z::Float64;
                        w0::Float64=-1.,
                        wprime::Float64 = 0.,
                        rhoDE::Float64 = 1.)

  ans = rhoDE*((1/(1+z))^(-(3+3*w0+6*wprime)))*exp(-6*wprime*(1-1/(1+z)))
  return ans
end
