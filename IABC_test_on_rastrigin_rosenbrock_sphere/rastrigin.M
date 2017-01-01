
function ObjVal = rastrigin(Chrom,switch1);

% Dimension of objective function
    Dim=size(Chrom,2);
   
% Compute population parameters
   [Nind,Nvar] = size(Chrom);

      % function 6, Dim*A + sum of (xi^2 - A*cos(Omega*xi)) for i = 1:Dim (Dim=20)
      % n = Dim, -5.12 <= xi <= 5.12
      % global minimum at (xi)=(0) ; fmin=0
      A = 10;
      Omega = 2 * pi;
      ObjVal = Dim * A + sum(((Chrom .* Chrom) - A * cos(Omega * Chrom))')';
      % ObjVal = diag(Chrom * Chrom');  % both lines produce the same
   % otherwise error, wrong format of Chrom
