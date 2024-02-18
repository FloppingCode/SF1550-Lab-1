%% 1c - Newtons metod
f = @(x) x.^2-8*x-12*sin(3*x+1)+19;
df = @(x) 2*x-8-36*cos(3*x+1);
starting_guesses = [1.97, 2.67, 3.9, 4.8, 6.2, 6.65]; 
tol = 1e-10;     
maxiter = 100; 

format long;

for guess = starting_guesses
    [root, iter, diffs] = newton(f, df, guess, tol, maxiter);
    fprintf('Nollställe: %12.8f\n', root);
    fprintf('Startgissning: %12.8f\n', guess);
    fprintf('Antal iterationer: %d\n', iter);
    fprintf('\n');
end

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
    
    % 1. Antalet korrekta siffror
    correct_digits = -log10(tol);
    actual_digits = -log10(abs(root - xold));
    fprintf('Antal korrekta siffror: %d\n', min(correct_digits, actual_digits));
end