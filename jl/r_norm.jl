function rfun(x, p)
    a = p[:a_r]; mur = p[:mu_r]; sigmar = p[:sigma_r]; mort = p[:m]
    a*exp(-(x-mur)^2/(2*sigmar^2)) - mort
end
function drdx(x, p)
  a = p[:a_r]; mur = p[:mu_r]; sigmar = p[:sigma_r]
  a*exp(-(x-mur)^2/(2*sigmar^2))*(-1/(sigmar^2))*(x-mur)
  # -((a * exp(-((-mur + x)^2/(2 * sigmar^2))) * (-mur + x))/sigmar^2)
end

function d2rdx(x, p)
  a = p[:a_r]; mur = p[:mu_r]; sigmar = p[:sigma_r]
  -(a/sigmar^2)*(exp(-(x-mur)^2/(2*sigmar^2))*(1-((x-mur)^2/(sigmar^2))))
end
