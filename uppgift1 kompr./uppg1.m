%% 1a - plotta f(x)

f = @(x) x.^2 - 8*x - 12*sin(3*x + 1) + 19;

figure(1);
fplot(f, [0, 5]);
ylim([-25, 25]);

xlabel('x');
ylabel('f(x)');

%% 1b - Fixpunktiteration
format long; 

x0_values = [1.97, 2.67, 3.9, 4.8, 6.2, 6.65];
tau = 1e-10;

for i = 1:length(x0_values)
    xit = fixpunkt(x0_values(i), tau);
    
    % Skriv ut de sista tio värdena på |xn+1 - xn| för den första roten
    if i == 1
        fprintf('De sista tio värdena på |xn+1 - xn|:\n');
        disp(abs(diff(xit(end-9:end))));
    end
    
    % Skriv ut nollstället och antalet iterationer utanför funktionen
    fprintf('För startgissning x0 = %f:\n', x0_values(i));
    if isempty(xit)
        fprintf('Nollstället kunde inte hittas.\n');
    else
        fprintf('Hittades nollstället %.20f med %d iterationer.\n', xit(end), length(xit));
    end
    
    % Skriv ut iterationer, approximationer och differenser för första roten
    if i == 1
        fprintf('Utskrift av iterationer, approximationer och differenser för första roten:\n');
        for j = 1:length(xit)
            fprintf('%4d      %12.8f   %12.8e\n', j, xit(j), abs(xit(j) - xit(max(1, j-1))));
        end
    end
    
    fprintf('--------------------------------------------\n');
end

%% 1c - Newton
format long;

x0_values = [1.97, 2.67, 3.9, 4.8, 6.2, 6.65];
tau = 1e-10;

for i = 1:length(x0_values)

    if i == 1
        fprintf('De sista tio värdena på |xn+1 - xn|:\n');
        disp(abs(diff(xit)));
    end

    xit = newton(x0_values(i), tau);
    
    fprintf('För startgissning x0 = %f:\n', x0_values(i));
    
    if isempty(xit)
        fprintf('Nollstället kunde inte hittas.\n');
    else
        fprintf('Hittades nollstället %f med antalet iterationer %d.\n', xit(end), length(xit));
    end
    
    fprintf('--------------------------------------------\n');
end


%% 1d - konvergensplottar
format long;

% Val av startgissning
x0 = 1.97;
tau = 1e-10;

% Referenslösning med Newtons metod
start_xref = newton(x0, 1e-15);
xref = start_xref(end);

xit_fixpunkt = fixpunkt(x0, tau);
xit_newton = newton(x0, tau);

% Subtrahera xref från båda vektorerna
diffs_newton = abs(xit_newton - xref);
diffs_fixpunkt = abs(xit_fixpunkt - xref);

% Plotta differenserna för både Newtons metod och fixpunktiteration
figure;
semilogy(1:length(diffs_newton), abs(diffs_newton), 'o-', 'DisplayName', 'Newton');
hold on;
semilogy(1:length(diffs_fixpunkt), abs(diffs_fixpunkt), 'x-', 'DisplayName', 'Fixpunkt');
xlabel('Iteration n');
ylabel('|xn - xref|');
title('Felet efter varje iteration för Newtons metod och fixpunktiteration');
legend('Location', 'Best');
grid on;

% Plotta konvergensordningen för både Newtons metod och fixpunktiteration
figure;
loglog(1:length(diffs_newton)-1, abs(diffs_newton(2:end)./diffs_newton(1:end-1)), 'o-', 'DisplayName', 'Newton');
hold on;
loglog(1:length(diffs_fixpunkt)-1, abs(diffs_fixpunkt(2:end)./diffs_fixpunkt(1:end-1)), 'x-', 'DisplayName', 'Fixpunkt');
xlabel('Iteration n');
ylabel('Konvergensordning');
title('Konvergensordning för Newtons metod och Fixpunktiteration');
legend('Location', 'Best');
grid on;
function xit = fixpunkt(x0, tau)
    g = @(x) 1/19 * (x^2 + 11*x - 12*sin(3*x+1)) + 1; 
    
    xold = x0;           
    diff = 1;           
    iter = 0;
    maxiter = 200;
    
    while diff > tau
        iter = iter + 1;           
        x = g(xold);               
        diff = abs(x - xold);      
        xold = x;                  
        
        xit(iter) = x;  
    end
    
    if iter == maxiter
        warning('Maximala antalet iterationer har uppnåtts%f.', xold);
        xit = []; 
    end
end

function xit = newton(x0, tau)
    f = @(x) x.^2 - 8*x - 12*sin(3*x + 1) + 19;
    df = @(x) 2*x - 8 - 36*cos(3*x+1);

    xold = x0;           
    diff = 1;           
    iter = 0;
    maxiter = 200;
        
    while diff > tau && iter < maxiter
        iter = iter + 1;           
        x = xold - (f(xold)/df(xold));               
        diff = abs(x - xold);  
        xold = x;                  
        
        xit(iter) = x;  
        
        %fprintf('%4d      %12.8f   %12.8e\n', iter, x, diff); %För att
        %visa att antalet korrkta siffror dubblas i varje iteration
    end
    
    fprintf('\n');
    
    if iter == maxiter
        warning('Maximala antalet iterationer har uppnåtts %f.', xold);
        xit = []; 
    end
end

