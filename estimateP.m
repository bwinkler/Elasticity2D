function pEst = estimateP( ds, z )
  pEst = ds.C \ (ds.B' * (ds.Q' * z));
end
