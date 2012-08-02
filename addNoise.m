function ud = addNoise( u, noisePercent, type )

  u1 = u(1:2:end);
  u2 = u(2:2:end);

  % unorms = zeros(size(u1));
  % for i = [1:size(u1)]
  %   unorms(i,1) = norm([u1(i), u2(i) ]);
  % end

  % a =  -noisePercent * max( unorms );
  % b = -a;
  % r = (a + (b-a) .* rand(size(u1)));
  % err_norm = norm(r);
  % ud1 = u1 + 0.5 * r;
  % ud2 = u2 + 0.5 * r;

  % ud = [];
  
  a = -noisePercent * max( u1 );
  b = -a;
  r = (a + (b-a) .* rand(size(u1)));
  ud1 = u1 + r;
  

  a = -noisePercent * max( u2);
  b = -a;
  r = (a + (b-a) .* rand(size(u2)));

  ud2 = u2 + r;

  ud = [];
  for i = [1:size(ud1)]
    ud = [ ud; ud1(i); ud2(i) ];
  end
end

  
