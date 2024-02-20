% Kod för uppgift 2...
%% 2a - Förberedande arbete

% Ekvationssystemet definieras enligt de givna ekvationerna
% |x1-a|^2 = ra^2 (1)
% |x2-b|^2 = rb^2 (2)
% (x1-x2)*(x1-a)=0 (3)
% (x1-x2)*(x2-b)=0 (4)

% För att skriva om systemet till F(x1, y1, x2, y2)=0 behöver vi definiera
% fyra vektorer F1 till F4 enligt följande:

% F1(x1, y1, x2, y2) = (x1-a(1))^2 + (y1-a(2))^2 - ra^2
% F2(x1, y1, x2, y2) = (x2-b(1))^2 + (y2-b(2))^2 - rb^2
% F3(x1, y1, x2, y2) = (x1-x2)*(x1-a(1))
% F4(x1, y1, x2, y2) = (x1-x2)*(x2-b(2))

% Matrisen F representerar vektorerna F1 till F4 för en given vektor x0 av variabler:
F = [((x0(1) - a(1))^2 + (x0(2) - a(2))^2 - ra^2);
     ((x0(3) - b(1))^2 + (x0(4) - b(2))^2 - rb^2);
     (x0(1) - x0(3))*(x0(1) - a(1)) + (x0(2) - x0(4))*(x0(2) - a(2));
     (x0(1) - x0(3))*(x0(3) - b(1)) + (x0(2) - x0(4))*(x0(4) - b(2))];

% Jakobimatrisen J blir då:
J = [2*(x0(1) - a(1)), 2*(x0(2) - a(2)), 0, 0;
     0, 0, 2*(x0(3) - b(1)), 2*(x0(4) - b(2));
     (x0(1) - a(1)) + (x0(1) - x0(3)), (x0(2) - a(2)) + (x0(2) - x0(4)), -(x0(1) - a(1)), -(x0(2) - a(2));
     -(x0(3) - b(1)), -(x0(4) - b(2)), (x0(1) - x0(3)) + (x0(3) - b(1)), (x0(2) - x0(4)) + (x0(4) - b(2))];

%% 2b - Löser ekvationssystemet (punkter)

x0 = [-1, 4, 1.5, 1.5] % Finns 4 lösningar, men om dessa startpunkter väljs visas 1 av dessa
ra = 1.5
rb = 0.8
a = [-1.5, 3.0]
b = [1.0, 1.0]

[Xrot, iter] = punkter(x0, ra, rb, a, b);

disp('Lösning:');
disp(['(x1, y1) = (' num2str(Xrot(1)) ', ' num2str(Xrot(2)) ')']);
disp(['(x2, y2) = (' num2str(Xrot(3)) ', ' num2str(Xrot(4)) ')']);

fprintf('Antalet iterationer: ')
disp(iter);

plot_2b(Xrot, a, b, ra, rb);


%% 2c - Beräknar snörets längd (langd)
format long;

a = [-1; 1.5];
b = [3.0; 0.5];
c = [0.0; -2.0];
ra = 1.0;
rb = 1.2;
rc = 1.7;
    
[L, X] = langd(ra, rb, rc, a, b, c);
disp('Snörets längd: ');
disp(L)
%disp(X)
plot_2c(L, X, a, b, c, ra, rb, rc);


%% 2d - Störningsanalys 
a = [-1; 1.5];
b = [3.0; 0.5];
c = [0.0; -2.0];
ra = 1.0;
rb = 1.2;
rc = 1.7;

osakerhet(ra, rb, rc, a, b, c)
%% Funktioner

%% Funktioner

function [L, X] = langd(ra, rb, rc, a, b, c)
    x0_AC = [-1.97; 1.28; -1.7; -2.1];
    x0_CB = [1.7; -2.2; 4.0; -0.2;];
    x0_BA = [3.4; 0.8; -0.8; 2.6];

    AC = punkter(x0_AC, ra, rc, a, c);
    CB = punkter(x0_CB, rc, rb, c, b);
    BA = punkter(x0_BA, rb, ra, b, a);

    %Beräkning av snörets längd
    length_AC = norm([AC(1) AC(2)] - [AC(3) AC(4)]);
    length_CB = norm([CB(1) CB(2)] - [CB(3) CB(4)]);
    length_BA = norm([BA(1) BA(2)] - [BA(3) BA(4)]);
    segment_A = circle_segment_length(a, [AC(1) AC(2)], [BA(3) BA(4)]);
    segment_B = circle_segment_length(b, [BA(1) BA(2)], [CB(3) CB(4)]);
    segment_C = circle_segment_length(c, [CB(1) CB(2)], [AC(3) AC(4)]);

    L = length_AC + length_BA + length_CB + segment_A + segment_B + segment_C;
    X = [BA(1:2), CB(1:2), AC(1:2); BA(3:4), CB(3:4), AC(3:4)];
end

function segment_length = circle_segment_length(center_coords, segment_coord1, segment_coord2)
    radius  = norm(segment_coord2 - center_coords');
    segment_distance = norm(segment_coord2 - segment_coord1);    
    theta = acos((-segment_distance^2 + 2*(radius^2))/(2*(radius^2)));
    segment_length = theta * radius;
end

function [Xrot, iter] = punkter(x0, ra, rb, a, b)

    % Startpunkter
    tol = 1e-10;

    % Newtons metod
    iteration = 0;
    while true
        F = [((x0(1) - a(1))^2 + (x0(2) - a(2))^2 - ra^2);
             ((x0(3) - b(1))^2 + (x0(4) - b(2))^2 - rb^2);
             (x0(1) - x0(3))*(x0(1) - a(1)) + (x0(2) - x0(4))*(x0(2) - a(2));
             (x0(1) - x0(3))*(x0(3) - b(1)) + (x0(2) - x0(4))*(x0(4) - b(2))];

        J = [2*(x0(1) - a(1)), 2*(x0(2) - a(2)), 0, 0;
             0, 0, 2*(x0(3) - b(1)), 2*(x0(4) - b(2));
             (x0(1) - a(1)) + (x0(1) - x0(3)), (x0(2) - a(2)) + (x0(2) - x0(4)), -(x0(1) - a(1)), -(x0(2) - a(2));
             -(x0(3) - b(1)), -(x0(4) - b(2)), (x0(1) - x0(3)) + (x0(3) - b(1)), (x0(2) - x0(4)) + (x0(4) - b(2))];

        x1 = x0 - J\F;

        if max(abs(x1 - x0)) < tol
            break;
        end

        x0 = x1;
        iteration = iteration + 1;
    end

    Xrot = x1;
    iter = iteration;

end

function osakerhet(ra, rb, rc, a, b, c)
    num_trials = 10;
    delta_values = 0.01 * (2 * rand(num_trials, 1) - 1); % Random delta values between -0.01 and 0.01
    
    [L, X] = langd(ra, rb, rc, a, b, c);
    params = {ra, rb, rc, a, b, c};

    for i = 1:numel(params)
        uncertainty_sum = 0;

        for j = 1:num_trials
            perturbed_params = params;
            
            % Viktigt att veta om man ska lägga till delta på skalär eller
            % vektor
            if isscalar(params{i})
                perturbed_params{i} = perturbed_params{i} + delta_values(j);
            else
                perturbed_params{i} = perturbed_params{i} + delta_values(j) * ones(size(params{i}));
            end

            [perturbed_L, ~] = langd(perturbed_params{1}, perturbed_params{2}, perturbed_params{3}, perturbed_params{4}, perturbed_params{5}, perturbed_params{6});

            uncertainty = (perturbed_L - L) / delta_values(j);
            uncertainty_sum = uncertainty_sum + uncertainty;
        end

        average_uncertainty = uncertainty_sum / num_trials;
        fprintf('Genomsnittlig osäkerhet för parameter %d: %f\n', i, average_uncertainty);
    end
    fprintf('\n')
end

function plot_2b(Xrot, a, b, ra, rb)

    % Plottning av cirklarna och snöret
    v = linspace(0, 2*pi, 100);
    circleA_x = a(1) + ra * cos(v);
    circleA_y = a(2) + ra * sin(v);

    circleB_x = b(1) + rb * cos(v);
    circleB_y = b(2) + rb * sin(v);

    line_x = [Xrot(1), Xrot(3)];
    line_y = [Xrot(2), Xrot(4)];


    figure;
    hold on;
    plot(circleA_x, circleA_y, 'r', 'LineWidth', 1);
    plot(circleB_x, circleB_y, 'b', 'LineWidth', 1);
    plot(line_x, line_y, 'g', 'LineWidth', 2);
    xlabel('x');
    ylabel('y');
    legend('Cirkel A', 'Cirkel B', 'Snöre', 'x1', 'x2');
    axis equal;
    grid on;
    hold off;
end

function plot_2c(L, X, a, b, c, ra, rb, rc)
    
    % Plottning av cirklarna och snöret
    v = linspace(0, 2*pi, 100);

    circleA_x = a(1) + ra * cos(v);
    circleA_y = a(2) + ra * sin(v);
    circleB_x = b(1) + rb * cos(v);
    circleB_y = b(2) + rb * sin(v);
    circleC_x = c(1) + rc * cos(v);
    circleC_y = c(2) + rc * sin(v);

    point_A = X(1:2, :);
    point_B = X(3:4, :);
    
    figure;
    hold on;
    
    plot([point_A(1, 1), point_B(1, 1)], [point_A(2, 1), point_B(2, 1)], 'g', 'LineWidth', 2);
    plot([point_A(1, 2), point_B(1, 2)], [point_A(2, 2), point_B(2, 2)], 'b', 'LineWidth', 2);
    plot([point_A(1, 3), point_B(1, 3)], [point_A(2, 3), point_B(2, 3)], 'r', 'LineWidth', 2);
    
    xlabel('x');
    ylabel('y');
    title('Snöre, tre cirklar');
    axis equal;
    grid on;
    hold off;

    hold on;
    axis equal;

    % Plottning av cirklarna och snöret
    plot(circleA_x, circleA_y, 'r', 'LineWidth', 1);
    plot(circleB_x, circleB_y, 'b', 'LineWidth', 1);
    plot(circleC_x, circleC_y, 'y', 'LineWidth', 1);
end


