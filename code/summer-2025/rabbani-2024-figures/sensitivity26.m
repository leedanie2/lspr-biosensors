% calculating S using formula 26 of Rabbani 2024 for 2 layer Au NPs
function sens26 = sensitivity26(lambda_star, l, t)
    n0 = 1.33; % NPs in suspension in water
    n = 1.40; % analyte layer RI
    neff = n + (n0 - n) * (l/(l+t)) ^ 3; % effective refractive index
    % fitting constants
    m = -0.0719;
    eps0 = 33.6;
   
    sens26 = (2/(m * neff)) * (m * lambda_star + eps0) * (1 - (l/(l+t)) ^ 3);