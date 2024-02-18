%% 1b - Fixpunktiteration

f = @(x) 1/19 * (x^2 + 11*x - 12*sin(3*x+1)) + 1;
it_f = @(x) 1/19 * (x^2 + 11*x - 12*sin(3*x+1)) + 1;

starting_guesses = [1.97];  
tol = 1e-10;     
maxiter = 100;   

for guess = starting_guesses
    fixpunkt_result = fixpunkt(it_f, guess, tol, maxiter);
    fprintf('Fixpunkt: ');
    disp(fixpunkt_result);
end

function x = fixpunkt(g, starting_guess, tol, maxiter)
    xold = starting_guess;  
    diff = 1;               
    iter = 0;               
    
    fprintf('Startgissning: %f\n', xold);
    
    while diff > tol && iter < maxiter
        iter = iter + 1;           
        x = g(xold);               
        diff = abs(x - xold);      
        xold = x;                  

        fprintf('%4d      %12.8f   %12.8e\n', iter, x, diff);
    end
    fprintf('\n');
    
    if iter == maxiter
        warning('Maximala antalet iterationer har uppnåtts utan att nollställets hittats för startgissning %f.', xold);
    end
end