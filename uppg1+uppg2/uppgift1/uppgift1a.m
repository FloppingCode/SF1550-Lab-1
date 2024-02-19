%% 1a - Plotta f(x)

f = @(x) x.^2 - 8*x - 12*sin(3*x + 1) + 19;

fplot(f, [0, 5]);
ylim([-25, 25]);

xlabel('x');
ylabel('f(x)');
