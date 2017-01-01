  function ObjVal = griewank(Chrom,switch1);
   [Nind, Nvar] = size(Chrom);
   nummer = repmat(1:Nvar, [Nind 1]);
   ObjVal = sum(((Chrom.^2) / 4000)')' - prod(cos(Chrom ./ sqrt(nummer))')' + 1; 
  