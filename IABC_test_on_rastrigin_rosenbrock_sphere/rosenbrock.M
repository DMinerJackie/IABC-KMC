function ObjVal = rosenbrock(Chrom,switc);

% Dimension of objective function

    Dim=size(Chrom,2);
   
% Compute population parameters
   [Nind,Nvar] = size(Chrom);


      % function 11, sum of 100* (x(i+1) -xi^2)^2+(1-xi)^2 for i = 1:Dim (Dim=10)
      % n = Dim, -10 <= xi <= 10
      % global minimum at (xi)=(1) ; fmin=0
      Mat1 = Chrom(:,1:Nvar-1);
      Mat2 = Chrom(:,2:Nvar);
     
      if Dim == 2
         ObjVal = 100*(Mat2-Mat1.^2).^2+(1-Mat1).^2;
      else
         ObjVal = sum((100*(Mat2-Mat1.^2).^2+(1-Mat1).^2)')';
      end   
  

% End of function

