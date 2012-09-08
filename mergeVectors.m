function u = mergeVectors(u1,u2)
  u = [];
  for i = [1:size(u1)]
    u = [ u; u1(i); u2(i) ];
  end
end  
