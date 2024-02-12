 % Startpunkter
a = [-1.5; 3.0];
b = [1.0; 1.0];
ra = 1.5;
rb = 0.8;

solution = punkter(a, b, ra, rb);

disp('Lösning:');
disp(['(x1, y1) = (' num2str(solution(1)) ', ' num2str(solution(2)) ')']);
disp(['(x2, y2) = (' num2str(solution(3)) ', ' num2str(solution(4)) ')']);


% Plottning av cirklarna och snöret
v = linspace(0, 2*pi, 100);
circleA_x = a(1) + ra * cos(v);
circleA_y = a(2) + ra * sin(v);

circleB_x = b(1) + rb * cos(v);
circleB_y = b(2) + rb * sin(v);

line_x = [solution(1), solution(3)];
line_y = [solution(2), solution(4)];

figure;
hold on;
plot(circleA_x, circleA_y, 'r', 'LineWidth', 1);
plot(circleB_x, circleB_y, 'b', 'LineWidth', 1);
plot(line_x, line_y, 'g', 'LineWidth', 2);
scatter(solution(1), solution(2), 'b', 'filled');
scatter(solution(3), solution(4), 'r', 'filled');
xlabel('x');
ylabel('y');
legend('Cirkel A', 'Cirkel B', 'Snöre', 'x1', 'x2');
axis equal;
grid on;
hold off;


function solution = punkter(a, b, ra, rb)
    x0 = [0; 0; 0; 0];

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

    solution = x1;
end
