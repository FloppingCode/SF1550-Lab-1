%% 1a - plotta f(x)

f = @(x) x.^2 - 8*x - 12*sin(3*x + 1) + 19;

figure(1);
fplot(f, [-7, 20]);
ylim([-25, 25]);

xlabel('x');
ylabel('f(x)');

%% 1b - Fixpunktiteration
format long; 

x0_values = [1.97, 2.67, 3.9, 4.8, 6.2, 6.65];
tau = 1e-10;

for i = 1:length(x0_values)
    xit = fixpunkt(x0_values(i), tau);
    
    fprintf('För startgissning x0 = %f:\n', x0_values(i));
    if isempty(xit)
        fprintf('Nollstället kunde inte hittas.\n');
    else
        fprintf('Nollstället %.20f hittades med %d iterationer.\n', xit(end), length(xit));
        if i == 1
            fprintf('De sista tio värdena på |xn+1 - xn|:\n');
            diffs = abs(diff(xit(end-9:end)));
            for k = 1:length(diffs)
                fprintf('%d. %.55f\n', k, diffs(k));
            end
        end
    end
    
    % Skriver ut iterationer, approximationer och differenser för första roten
    if i == 1
        %fprintf('Utskrift av iterationer, approximationer och differenser för första roten:\n');
        for j = 1:length(xit)
            %fprintf('%4d      %12.8f   %12.8e\n', j, xit(j), abs(xit(j) - xit(max(1, j-1))));
        end
    end
    
    fprintf('--------------------------------------------\n');
end

%% 1c - Newton

format long;

x0_values = [1.97, 2.67, 3.9, 4.8, 6.2, 6.65];
tau = 1e-10;

for i = 1:length(x0_values)
    xit = newton(x0_values(i), tau);

    fprintf('För startgissning x0 = %f:\n', x0_values(i));

    if isempty(xit)
        fprintf('Nollstället kunde inte hittas.\n');
    else
        fprintf('Hittades nollstället %f med %d iterationer.\n', xit(end), length(xit));
        
        if i == 1
            diffs = abs(diff(xit));
            fprintf('Värdena på |xn+1 - xn| för den första roten:\n');
            for j = 1:length(diffs)
                fprintf('%.55f\n', diffs(j));
            end
        end
    end
    
    fprintf('--------------------------------------------\n');
end



%% 1d - konvergensplottar

format long;

x0 = 2;
tau = 1e-10;

% Referenslösning med Newtons metod
start_xref = newton(x0, 1e-15);
xref = start_xref(end);

xit_fixpunkt = fixpunkt(x0, tau);
xit_newton = newton(x0, tau);

diffs_newton = abs(xit_newton - xref);
diffs_fixpunkt = abs(xit_fixpunkt - xref);

% Plotta differenserna för både Newtons metod och fixpunktiteration
figure;
semilogy(1:length(diffs_newton), abs(diffs_newton), 'o-', 'DisplayName', 'Newton');
hold on;
semilogy(1:length(diffs_fixpunkt), abs(diffs_fixpunkt), 'x-', 'DisplayName', 'Fixpunkt');
xlabel('Iteration n');
ylabel('|xn - xref|');
title('Felet efter varje iteration');
legend('Location', 'Best');
grid on;

xit_fixpunkt = fixpunkt(x0, tau);
xit_newton = newton(x0, tau);

% Beräkna felet |en+1| och |en|och plottar
errors_fixpunkt = abs(diff(xit_fixpunkt));
errors_newton = abs(diff(xit_newton));
figure;
loglog(errors_fixpunkt(1:end-1), errors_fixpunkt(2:end), 'o-', 'DisplayName', 'Fixpunkt (p = 1)');
hold on;
loglog(errors_newton(1:end-1), errors_newton(2:end), 'x-', 'DisplayName', 'Newton (p = 2)');

% Plotta linjer med exakta lutningarna 1 och 2
loglog(errors_fixpunkt(1:end-1), errors_fixpunkt(1:end-1), '--', 'DisplayName', 'Exakt lutning (p = 1)');
loglog(errors_newton(1:end-1), errors_newton(1:end-1).^2, '--', 'DisplayName', 'Exakt lutning (p = 2)');

xlabel('|en|');
ylabel('|en+1|');
title('|en+1| som funktion av |en|');
legend('Location', 'Best');
grid on;

%% Funktioner
%% Funktioner
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


