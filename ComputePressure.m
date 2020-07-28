function [P] = ComputePressure(rho,c,T_TR,T_RT,T_RM,T_TM,t1,t2,f1,f2,d1,d2,nT,nR,nM,type1,type2)

  % 1 IS UP
  % 2 IS DOWN
  PT0 = zeros(nM,1);
  PT1 = zeros(nM,1);
  PT2 = zeros(nM,1);
  PT3 = zeros(nM,1);
  PT4 = zeros(nM,1);
  PT5 = zeros(nM,1);
  PT6 = zeros(nM,1);

  PR0 = zeros(nM,1);
  PR1 = zeros(nM,1);

  if f1 == 0
    C1 = 0;
    A1 = 0;
    U1 = zeros(nT,1);
  else
    wL1 = c/f1;
    A1 = (1j./wL1);
    omega1 = 2 * pi * f1;
    C1 = omega1*rho*c/wL1;
    U1 = ones(nT ,1).*d1.*exp(-1j*(omega1*t1));
  end

  if f2 == 0
    C2 = 0;
    A2 = 0;
    U2 = zeros(nR,1);
  else
    wL2 = c/f2;
    A2 = (-1j./wL2);
    omega2 = 2 * pi * f2;
    C2 = omega2*rho*c/wL2;
    U2 = ones(nR ,1).*d2.*exp(-1j*(omega2*t2));
  end

  if (strcmp(type1, "Transducer") == 1 && strcmp(type2, "Reflector") == 1)
    PT0 = (C1)*T_TM*U1;
    PT1 = (C1)*(A1)*T_RM*T_TR*U1;
    PT2 = (C1)*(A1^2)*T_TM*T_RT*T_TR*U1;
    PT3 = (C1)*(A1^3)*T_RM*T_TR*T_RT*T_TR*U1;
    PT4 = (C1)*(A1^4)*T_TM*T_RT*T_TR*T_RT*T_TR*U1;
    PT5 = (C1)*(A1^5)*T_RM*T_TR*T_RT*T_TR*T_RT*T_TR*U1;
    PT6 = (C1)*(A1^6)*T_TM*T_RT*T_TR*T_RT*T_TR*T_RT*T_TR*U1;

  elseif (strcmp(type1, "Array") == 1 && strcmp(type2, "Array") == 1)
    PT0 = (C1)*T_TM*U1;
    PR0 = (C2)*T_RM*U2;

    PT = [PT0];
    PR = [PR0];

  elseif (strcmp(type1, "Reflector") == 1 && strcmp(type2, "Array") == 1)
    PR0 = (C2)*T_RM*U2;
    PR1 = (C2)*(A2)*T_TM*T_RT*U2;

    PT = [PT0, PT1];
    PR = [PR0, PR1];

  elseif (strcmp(type1, "Array") == 1 && strcmp(type2, "Reflector") == 1)
    PT0 = (C1)*T_TM*U1;
    PT1 = (C1)*(A1)*T_RM*T_TR*U1;

    PT = [PT0, PT1];
    PR = [PR0, PR1];
  end

  % p2p_PT = peak2peak(real(PT))
  % p2p_PR = peak2peak(real(PR))
  %
  % figure(4)
  % plot(p2p_PT, '-r', 'LineWidth', 2)
  % figure(5)
  % plot(p2p_PR, '-b', 'LineWidth', 2)

  PT = [PT0 PT1 PT2 PT3 PT4 PT5 PT6];
  PR = [PR0 PR1];
  P = sum(PT,2) + sum(PR,2);

end
