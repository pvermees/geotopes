function se = S2se(S)
  covmat = (inv(S.R)*inv(S.R'))*S.normr^2/S.df;
  se = sqrt(diag(covmat));
end