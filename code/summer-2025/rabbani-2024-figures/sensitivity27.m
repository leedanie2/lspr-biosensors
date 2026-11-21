% calculating S using formula 27 of Rabbani 2024 for 2 layer Au NPs
function sens27 = sensitivity27(lambda0, l, t)
    n0 = 1.33; % NPs in suspension in water
    n = 1.40; % analyte layer RI
    neff = n + (n0 - n) * (l/(l+t)) ^ 3; % effective refractive index
    % fitting constants
    m = -0.0719;
    eps0 = 33.6;

    sens27 = (2 * neff/(m * n0 ^ 2)) * (m * lambda0 + eps0) * (1 - (l/(l+t)) ^ 3);