% Define the range of X and Y values for the phase plane
X_range = linspace(0, 4, 100);
Y_range = linspace(0, 4, 100);

% Create a grid of X and Y values
[X, Y] = meshgrid(X_range, Y_range);

% Calculate the derivatives at each point in the grid
n = 100;
dXdt = 0.2 + (3./(1 + Y.^n)) - X;
dYdt = 0.2 + (3./(1 + X.^n)) - Y;

% Normalize the derivatives for plotting
magnitude = sqrt(dXdt.^2 + dYdt.^2);
dXdt = dXdt ./ magnitude;
dYdt = dYdt ./ magnitude;

% Create a phase plane plot
figure;
quiver(X, Y, dXdt, dYdt, 0.5);
xlabel('X');
ylabel('Y');
title('Mutually Inhibiting Circuit Phase Plane');

% Find equilibrium points
[X_eq, Y_eq] = find(magnitude < 0.1);

% Plot the equilibrium points
hold on;
scatter(X_range(X_eq), Y_range(Y_eq), 'r', 'filled');

% Add contour lines for dX/dt = 0 and dY/dt = 0
contour(X, Y, dXdt, [0 0], 'b');  % Contour for dX/dt = 0
contour(X, Y, dYdt, [0 0], 'g');  % Contour for dY/dt = 0

legend('Phase portrait lines','Equilibrium Points', 'dY/dt = 0', 'dX/dt = 0');
