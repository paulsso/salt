function T_TM = MakeTTM(type1, type2, rho, c, f1, f2, r_nm)
% -------------------------------------------------------------------------
% % Computation transfer matrices
% -------------------------------------------------------------------------
  Si = 1e-6;
  Sn = 1e-6;

  if f1 == 0
    kk1 = 0;
  else
    wL1 = c/f1;
    kk1 = 2 * pi / wL1;
  end

  if f2 == 0
    kk2 = 0;
  else
    wL2 = c/f2;
    kk2 = 2 * pi / wL2;
  end

  if strcmp(type1, "Array") == 1 && strcmp(type2, "Array") == 1
      T_TM = Si*exp(-1j*kk2*r_nm)./r_nm; % Transfer matrix down/medium
  elseif strcmp(type1, "Array") == 1 && strcmp(type2, "Reflector") == 1
      T_TM = Si*exp(-1j*kk1*r_nm)./r_nm; % Transfer matrix down/medium
  elseif strcmp(type1, "Reflector") == 1 && strcmp(type2, "Array") == 1
      T_TM = Si*exp(-1j*kk2*r_nm)./r_nm; % Transfer matrix down/medium
  elseif strcmp(type1, "Reflector") == 1 && strcmp(type2, "Reflector") == 1
      T_TM = 0*Si*exp(-1j*kk2*r_nm)./r_nm; % Transfer matrix down/medium
  end
end
