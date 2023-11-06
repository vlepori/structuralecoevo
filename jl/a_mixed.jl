# Asymm

function a0(x1,x2, p)
  c = p[:c]; v = p[:v]; k = p[:k]
  c * (1 - 1/(v* exp(-k * (x1 - x2)) + 1))
end
function a1(x1,x2, p)
  c = p[:c]; v = p[:v]; k = p[:k]
  -c* k* v*(exp(-k*(x1-x2)))*(1/((1+v*exp(-k*(x1-x2)))^2))
end

function a2(x1,x2, p)
  c = p[:c]; v = p[:v]; k = p[:k]
  (c * k^2 * v * exp(k *(x1 - x2)) * (exp(k * (x1 - x2)) - v))/(exp(k * (x1 - x2)) + v)^3
end

# Symm

function s0(x1,x2, p)
  a, sigma = p[:a], p[:sigma]
  a*exp(-((x1-x2)^2/(2*sigma^2)))
  # (2*pi)^(.5)*a*a*((sigma*sigma)/sqrt(sigma^2+sigma^2))*exp(-((x1-x2)^2/(2*(sigma^2+sigma^2))))
end
function s1(x1, x2, p)
  a, sigma = p[:a], p[:sigma]
  -(a*exp(-((x1-x2)^2/(2*sigma^2)))*(x1-x2))/(sigma^2)
  # (2*pi)^(.5)*a*a*((sigma*sigma)/sqrt(sigma^2+sigma^2))*exp(-((x1-x2)^2/(2*(sigma^2+sigma^2))))*-(1/(2*sigma^2))*(x1-x2)
end

function s2(x1,x2, p)
  a, sigma = p[:a], p[:sigma]
  (a*exp(-((x1-x2)^2/(2*sigma^2)))*(x1^2 - 2*x1*x2 + x2^2 - sigma^2))/(sigma^4)
  # -sqrt(2*pi)*a*a*((sigma*sigma)/(sigma^2+sigma^2)^(3/2))*(exp(-(x1-x2)^2/(4*sigma^2))*(1 - (x1-x2)^2/(2*sigma^2)))
end

# mixed

function alphafun(x1,x2, p)
  beta = p[:beta]
  beta*a0(x1,x2,p) + (1-beta)*s0(x1,x2,p)
end
function dadx(x1,x2, p)
  beta = p[:beta]
  beta*a1(x1,x2,p) + (1-beta)*s1(x1,x2,p)
end
function d2adx(x1,x2, p)
  beta = p[:beta]
  beta*a2(x1,x2,p) + (1-beta)*s2(x1,x2,p)
end
