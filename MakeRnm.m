function r_nm = MakeRnm(Vx, Vy, Vz, Mx, My, Mz)

  nT = length(Vx);
  nM = length(Mx);

  Ax = repmat(Vx',nM,1);
  Bx = repmat(Mx,1,nT);

  Ay = repmat(Vy',nM,1);
  By = repmat(My,1,nT);

  Az = repmat(Vz',nM,1);
  Bz = repmat(Mz,1,nT);

  A = Bx-Ax;
  B = By-Ay;
  C = Bz-Az;
  r_nm= sqrt(A.^2 + B.^2 + C.^2);
end
