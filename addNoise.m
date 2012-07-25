function [ud1, ud2, err_norm] = addNoise( u1, u2, noisePercent, type )

  unorms = zeros(size(u1));
  for i = [1:size(u1)]
    unorms(i,1) = norm([u1(i), u2(i) ]);
  end

  a =  -noisePercent * max( unorms );
  b = -a;
  r = (a + (b-a) .* rand(size(u1)));
  err_norm = norm(r);
  ud1 = u1 + 0.5 * r;
  ud2 = u2 + 0.5 * r;
end

  
