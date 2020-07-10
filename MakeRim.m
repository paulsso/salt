function r_im = MakeRim(Ux, Uy, Uz, Mx, My, Mz)

  nM = length(Mx);
  nR = length(Ux);

  Ax = repmat(Ux',nM,1);
  Bx = repmat(Mx,1,nR);

  Ay = repmat(Uy',nM,1);
  By = repmat(My,1,nR);

  Az = repmat(Uz',nM,1);
  Bz = repmat(Mz,1,nR);

  A = Bx-Ax;
  B = By-Ay;
  C = Bz-Az;
  r_im = sqrt(A.^2 + B.^2 + C.^2);
end
