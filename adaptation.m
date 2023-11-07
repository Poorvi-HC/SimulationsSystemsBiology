% Parameters
beta1 = 0.1;
beta2 = 0.1;
alpha1 = 0.1;
alpha2 = 0.1;
t = 0:1500;

% Define the time points when X is stepped up
step_x_time = t(21:125:end);  % Add more step times if needed

% Initial conditions
initial_state = [1.0; 1.0; 1.0];

% Solve the ODEs using 4th-order Runge-Kutta
solution = RungeKutta(t, initial_state, @model, beta1, beta2, alpha1, alpha2, step_x_time);

% Extract the variables
X = solution(:, 1);
Y = solution(:, 2);
Z = solution(:, 3);

% Plot the results
figure;

subplot(3, 1, 1);
plot(t, X);
xlabel('Time');
ylabel('[X]');
title('[X] vs time');

subplot(3, 1, 2);
semilogy(t, Y);  % Use semilogy for a logarithmic Y-axis
xlabel('Time');
ylabel('ln([Y])');
title('[Y] vs time');

subplot(3, 1, 3);
plot(t, Z);
xlabel('Time');
ylabel('[Z]');
title('[Z] vs time');

% Define the model function
function dy = model(state_vector, t, beta1, beta2, alpha1, alpha2, step_x_time)
    X = state_vector(1);
    Y = state_vector(2);
    Z = state_vector(3);

    % Calculate the derivatives
    dX = 0;
    if ismember(t, step_x_time)
        dX = 12;
    end

    dZ = (X * beta2 / Y) - Z * alpha2;
    dY = X * beta1 - Y * alpha1;

    dy = [dX; dY; dZ];
end

% Define the function for 4th-order Runge-Kutta
function y = RungeKutta(t, y0, f, beta1, beta2, alpha1, alpha2, step_x_time)
    n = length(t);
    y = zeros(n, length(y0));
    y(1, :) = y0;
    for i = 1:n-1
        h = t(i+1) - t(i);
        k1 = f(y(i, :), t(i), beta1, beta2, alpha1, alpha2, step_x_time).';  % Transpose k1
        k2 = f(y(i, :) + h/2 * k1, t(i) + h/2, beta1, beta2, alpha1, alpha2, step_x_time).';  % Transpose k2
        k3 = f(y(i, :) + h/2 * k2, t(i) + h/2, beta1, beta2, alpha1, alpha2, step_x_time).';  % Transpose k3
        k4 = f(y(i, :) + h * k3, t(i+1), beta1, beta2, alpha1, alpha2, step_x_time).';  % Transpose k4
        y(i+1, :) = (h/6) * (k1 + 2*k2 + 2*k3 + k4) + y(i, :);
    end
end