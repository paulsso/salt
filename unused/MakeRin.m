function r_in = MakeRin(Vx, Vy, Vz, Ux, Uy, Uz)

  nT = length(Vx);
  nR = length(Ux);

  Ax = repmat(Ux,1,nT);
  Bx = repmat(Vx',nR,1);

  Ay = repmat(Uy,1,nT);
  By = repmat(Vy',nR,1);

  Az = repmat(Uz,1,nT);
  Bz = repmat(Vz',nR,1);

  A = Bx-Ax;
  B = By-Ay;
  C = Bz-Az;
  r_in = sqrt(A.^2 + B.^2 + C.^2);

end
