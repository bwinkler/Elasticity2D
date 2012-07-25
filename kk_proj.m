function P = kk_proj(kku, kkl)
% KK_PROJ  Returns a projector onto the active set
% KK_PROJ(x, kku, kkl) projects x onto the active set with upper bound kku and
% lower bound kkl.
  if( max(kkl > kku) == 1  )
    error('Lower bound exceeds upper bound.');
  elseif( max(kku < kkl) == 1)
    error('Upper bound less than lower bound.');
  end

  P = @(x) max(kkl, min(kku,x));
