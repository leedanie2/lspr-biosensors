% Martinsson Figure 2(a)/2(e) Ag Nanosphere Sensitivity in solution vs substrate

n = zeros(1,10);
deltaLambdaSub = zeros(1,10);
deltaLambdaSol = zeros(1,10);

index = 1;
for i = 1.33:0.01:1.42
    n(index) = i;
    [ deltaLambdaSub(index), deltaLambdaSol(index) ] = sphereLambdaDrop(i, 400);
    index = index + 1;
end

%% 

plot(n, deltaLambdaSub, 'r'); hold on
plot(n, deltaLambdaSol, 'k'); hold off

xlim([1.33 1.42]);
ylim([0 30]);
xlabel( 'Refractive Index' );
ylabel( 'Δλmax / nm' );