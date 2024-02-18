%% 1d - Jämförelse av konvergens
f = @(x) x.^2-8*x-12*sin(3*x+1)+19;
df = @(x) 2*x-8-36*cos(3*x+1);
it_f = @(x) 1/19 * (x^2 + 11*x - 12*sin(3*x+1)) + 1;

% Referenslösning med Newtons metod
[reference_solution, ~, ~] = newton(f, df, 1.97, 1e-15, 100);

% Jämförelse av konvergensen
figure;
[x_newton, iter_newton, diffs_newton] = newton(f, df, 1.97, 1e-10, 100);
[root_fixpoint, iter_fixpoint, diffs_fixpoint] = fixpunkt(it_f, 1.97, 1e-10, 100);

% Plotta skillnaderna
plot(1:iter_newton, diffs_newton, '-o', 'DisplayName', 'Newton');
hold on;
plot(1:iter_fixpoint, diffs_fixpoint, '-o', 'DisplayName', 'Fixpunktiteration');
hold off;

title('Jämförelse av konvergens');
xlabel('Iteration');
ylabel('|xn - x*|');
legend('Newton', 'Fixpunktiteration');

function [root, iter, diffs] = newton(f, df, starting_guess, tol, maxiter)
    xold = starting_guess;  
    diff = 1;               
    iter = 0;               
    diffs = [];
    
    fprintf('Startgissning: %12.8f\n', xold);
    
    while diff > tol && iter < maxiter
        iter = iter + 1;     
        x = xold - (f(xold)/df(xold));               
        diff = abs(x - xold);  
        xold = x;                  
        diffs(iter) = diff;

        fprintf('%4d      %12.8f   %12.8e\n', iter, x, diff);
    end
    
    fprintf('\n');
    
    if iter == maxiter
        warning('Maximala antalet iterationer har uppnåtts %f.', xold);
    end
    
    root = x;  
end

function [root, iter, diffs] = fixpunkt(it_f, starting_guess, tol, maxiter)
    xold = starting_guess;  
    diff = 1;               
    iter = 0;               
    diffs = [];
    
    fprintf('Startgissning: %f\n', xold);
    
    while diff > tol && iter < maxiter
        iter = iter + 1;           
        x = it_f(xold);  % Använd it_f(xold) för fixpunktiterationen              
        diff = abs(x - xold);      
        xold = x;                  
        diffs(iter) = diff;

        fprintf('%4d      %12.8f   %12.8e\n', iter, x, diff);
    end
    fprintf('\n');
    
    if iter == maxiter
        warning('Maximala antalet iterationer har uppnåtts %f.', xold);
    end
    
    root = x;
end