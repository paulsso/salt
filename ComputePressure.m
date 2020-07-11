function [P] = ComputePressure(rho,c,T_TR,T_RT,T_RM,T_TM,t1,t2,f1,f2,d1,d2,nT,nR,depth1,depth2)

  % 1 IS UP
  % 2 IS DOWN
  T_TR = sparse(T_TR);
  T_RT = sparse(T_RT);
  T_RM = sparse(T_RM);
  T_TM = sparse(T_TM);
% TODO : ADD PHASE-CONFIGURATION
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
    A2 = (1j./wL2);
    omega2 = 2 * pi * f2;
    C2 = omega2*rho*c/wL2;
    U2 = -ones(nR ,1).*d2.*exp(-1j*(omega2*t2));
  end

  PT = (C1)*T_TM*U1...
    + (C1)*(A1)*T_RM*T_TR*U1...
    + (C1)*(A1^2)*T_TM*T_RT*T_TR*U1 ...
    + (C1)*(A1^3)*T_RM*T_TR*T_RT*T_TR*U1 ...
    + (C1)*(A1^4)*T_TM*T_RT*T_TR*T_RT*T_TR*U1...
    + (C1)*(A1^5)*T_RM*T_TR*T_RT*T_TR*T_RT*T_TR*U1...
    + (C1)*(A1^6)*T_TM*T_RT*T_TR*T_RT*T_TR*T_RT*T_TR*U1;
    
  PR = (C2)*T_RM*U2
    + (C2)*(A2)*T_TM*T_RT*U2...
    + (C2)*(A2^2)*T_RM*T_TR*T_RT*U2...
    + (C2)*(A2^3)*T_TM*T_RT*T_TR*T_RT*U2...
    + (C2)*(A2^4)*T_RM*T_TR*T_RT*T_TR*T_RT*U2...
    + (C2)*(A2^5)*T_TM*T_RT*T_TR*T_RT*T_TR*T_RT*U2...
    + (C2)*(A2^6)*T_RM*T_TR*T_RT*T_TR*T_RT*T_TR*T_RT*U2;

  P = PT + PR;
end
