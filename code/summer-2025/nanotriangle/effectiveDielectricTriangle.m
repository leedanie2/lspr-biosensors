function epseff = effectiveDielectricTriangle( nSub, nSolv, length, height)
    xA = (sqrt(3) * length / 4)/(sqrt(3) * length / 2 + 3 * height);
    epsSub = nSub ^ 2;
    epsSolv = nSolv ^ 2;
    epseff = xA * epsSub + (1 - xA) * epsSolv;
end