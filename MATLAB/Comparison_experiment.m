clear;
clc;
close all;

lightpink = [1, 0.59, 0.59];
lightblue = [0.62, 0.89, 1];


% parameters

beta = 1/4114;  %inverse temperature (fN times nm)
l_P = 33.9;     %persistence length (nm)
L = 103.7;      %DNA length (nm)


% Two-state model

ext = @(f) L.*((4/3) - 4 ./ (3 .* sqrt(f .* l_P .* beta + 1)) - ...
    (10 .* exp((900 ./ (beta .* f .* l_P)).^(1./4))) ./ ...
    (sqrt(f .* l_P .* beta) .* (exp((900 ./ (beta .* f .* l_P)).^(1./4)) - 1).^2) + ...
    ((f .* l_P .* beta).^1.62) ./ (3.55 + 3.8 .* (f .* l_P .* beta).^2.2));  %force-extension

betaF = @(x) 1./l_P.*(-L.^2./(4.*(x - L)) - x./(4) + x.^2./(2.*L));     %free energy F without force evaluated in x
betaG = @(f) betaF(ext(f)) - beta.* ext(f) .* f;                        %free energy G with force evaluated in f
betaG5 = @(f) betaF(5) - beta.* 5.* f;                                  %free energy G evaluated at f(x = 5 nm, which is the size of the protein)

DeltabetaG = @(f) -betaG(f) + betaG5(f);    %energy difference

tau = @(f) exp(DeltabetaG(f));      %looping time as predicted by two-state model
tau = @(f) tau(f)./tau(60);         %normalised at 60 fN

fplot(tau, 'LineWidth', 3, 'Color', lightblue, 'LineStyle', ':');

% Barrier escape model

f_range = [0:10:210];       %force range in fN
Lcut = 0;
x_c = 5;                    %protein size 5 nm

% Define the limits of integration
x_lower = x_c;
x_upper = L-Lcut;

% Define the function for the upper limit of y in terms of x
y_upper_func = @(x) x;

% Define the function for the upper limit of z in terms of y
z_upper_func = @(y) L-Lcut;

% Define the limits of integration for y and z as functions of x
y_lower_func = @(x) x_c;
z_lower_func = @(y) y;

% Define the numerical integration method
num_points = 100; % Number of points for numerical integration

% Initialize the result
result = 0;

tau = zeros(length(f_range),1);

for f_index = 1:length(f_range)
    result = 0;

    f = f_range(f_index);

    F = @(x) -log(x.^2.*(1 - (x./L).^2).^(-9./2).*exp(-3.*L./(4.*l_P.*(1 - (x./L).^2)))) - beta.*f.*x;


    % Define the function to integrate over x, y, and z
    fun = @(x, y, z) exp(-(F(x) - F(y) + F(z)));

    integrand_Z = @(x) exp(-F(x));
    Z = integral(integrand_Z, x_c, L-Lcut);


    % Perform the triple integral using numerical quadrature
    for i = 1:num_points
        x = x_lower + (x_upper - x_lower) * (i - 0.5) / num_points;

        y_lower = y_lower_func(x);
        y_upper = y_upper_func(x);

        for j = 1:num_points
            y = y_lower + (y_upper - y_lower) * (j - 0.5) / num_points;

            z_lower = z_lower_func(y);
            z_upper = z_upper_func(y);

            for k = 1:num_points
                z = z_lower + (z_upper - z_lower) * (k - 0.5) / num_points;

                % Evaluate the function and multiply by the volume element
                result = result + fun(x, y, z) * ((x_upper - x_lower) / num_points) * ((y_upper - y_lower) / num_points) * ((z_upper - z_lower) / num_points);
            end
        end
    end

    tau(f_index) = result./(Z);
    f

end

hold on
plot(f_range, tau./tau(7),'LineWidth',3,'LineStyle','--','Color',lightpink)

% Experimental data from Physical Review Letters 104, 258103 (2010),
% extracted using PlotDigitalizer
xx = [59.97, 77.81, 122.94, 152.91, 175.03, 183.06];
yy = [3.77, 5.15, 11.30, 14.72, 20.23, 30.98];

plot(xx, yy./yy(1), 'o', 'LineWidth', 2, 'LineStyle', 'none', 'MarkerSize', 8, 'Color', 'k');

set(gca, 'YScale', 'log');

lgnd = legend({'Two-state', 'Barrier escape', 'Experiment'},'Location','southeast','Interpreter','latex','FontSize',14);

box on
ax = gca;
ax.FontSize = 14; 

xlabel('$f$ [fN]', 'Interpreter','latex', 'FontSize',18)
ylabel('$\tau(f)/\tau(f=60$ fN$)$', 'Interpreter','latex', 'FontSize',18)

set(lgnd,'color','none');

xlim([25 215])

annotation('textbox',...
    [0.130768514244485 0.753195690889718 0.278703834758237 0.137756690062671],...
    'VerticalAlignment','cap',...
    'String',{'$L = 103.7$ nm','$l_P = 33.9$ nm'},...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');
