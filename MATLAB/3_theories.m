clear;
clc;
close all;

hfig = figure(1);

purple = [160/256, 156/256, 252/256];
darkbrown = [0.619607843137255 0.192156862745098 0.192156862745098];
lightbrown = [1, 0.59, 0.59];
lightblue = [0.62, 0.89, 1];

beta = 1;   %inverse temperature dimensionless units
l_P = 5;    %persistence length dimensionless units

%% READ SIMULATION DATA

% N = 11 (L/l_P = 2.2)
time0 = readtable("N11_f0.txt");
time0 = table2array(time0);
time0 = 0.005.*(time0);
t0 = mean(time0);
time0_1 = readtable("N11_f0_1.txt");
time0_1 = table2array(time0_1);
time0_1 = 0.005.*(time0_1);
t0_1 = mean(time0_1);
time0_2 = readtable("N11_f0_2.txt");
time0_2 = table2array(time0_2);
time0_2 = 0.005.*(time0_2);
t0_2 = mean(time0_2);
time0_3 = readtable("N11_f0_3.txt");
time0_3 = table2array(time0_3);
time0_3 = 0.005.*(time0_3);
t0_3 = mean(time0_3);
time0_4 = readtable("N11_f0_4.txt");
time0_4 = table2array(time0_4);
time0_4 = 0.005.*(time0_4);
t0_4 = mean(time0_4);
time0_5 = readtable("N11_f0_5.txt");
time0_5 = table2array(time0_5);
time0_5 = 0.005.*(time0_5);
t0_5 = mean(time0_5);
time0_6 = readtable("N11_f0_6.txt");
time0_6 = table2array(time0_6);
time0_6 = 0.005.*(time0_6);
t0_6 = mean(time0_6);
time0_7 = readtable("N11_f0_7.txt");
time0_7 = table2array(time0_7);
time0_7 = 0.005.*(time0_7);
t0_7 = mean(time0_7);
time0_8 = readtable("N11_f0_8.txt");
time0_8 = table2array(time0_8);
time0_8 = 0.005.*(time0_8);
t0_8 = mean(time0_8);
time0_9 = readtable("N11_f0_9.txt");
time0_9 = table2array(time0_9);
time0_9 = 0.005.*(time0_9);
t0_9 = mean(time0_9);
time1 = readtable("N11_f1.txt");
time1 = table2array(time1);
time1 = 0.005.*(time1);
t1 = mean(time1);

times11 = [time0, time0_1, time0_2, time0_3, time0_4, time0_5, time0_6, time0_7, time0_8, time0_9, time1];      %loop formation times all seeds
TT11 = [t0, t0_1, t0_2, t0_3, t0_4, t0_5, t0_6, t0_7, t0_8, t0_9, t1];                                          %loop formation time averages
t0_on11 = t0;                                                                                                   %loop formation time average at f = 0
TT11 = TT11.*1;                                                                                                 %C = 1 (constant prefactor)
FF11 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];                                                     %forces


% N = 23 (L/l_P = 4.6)
time0 = readtable("N23_f0.txt");
time0 = table2array(time0);
time0 = 0.005.*(time0);
t0 = mean(time0);
time0_1 = readtable("N23_f0_1.txt");
time0_1 = table2array(time0_1);
time0_1 = 0.005.*(time0_1);
t0_1 = mean(time0_1);
time0_2 = readtable("N23_f0_2.txt");
time0_2 = table2array(time0_2);
time0_2 = 0.005.*(time0_2);
t0_2 = mean(time0_2);
time0_3 = readtable("N23_f0_3.txt");
time0_3 = table2array(time0_3);
time0_3 = 0.005.*(time0_3);
t0_3 = mean(time0_3);
time0_4 = readtable("N23_f0_4.txt");
time0_4 = table2array(time0_4);
time0_4 = 0.005.*(time0_4);
t0_4 = mean(time0_4);
time0_5 = readtable("N23_f0_5.txt");
time0_5 = table2array(time0_5);
time0_5 = 0.005.*(time0_5);
t0_5 = mean(time0_5);
time0_6 = readtable("N23_f0_6.txt");
time0_6 = table2array(time0_6);
time0_6 = 0.005.*(time0_6);
t0_6 = mean(time0_6);
time0_7 = readtable("N23_f0_7.txt");
time0_7 = table2array(time0_7);
time0_7 = 0.005.*(time0_7);
t0_7 = mean(time0_7);


times23 = [time0, time0_1, time0_2, time0_3, time0_4, time0_5, time0_6];
TT23 = [t0, t0_1, t0_2, t0_3, t0_4, t0_5, t0_6]; 
t0_on23 = t0;
TT23 = TT23.*2;                                                                            %here C = 2
FF23 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6];


% N = 30 (L/l_P = 6)
time0 = readtable("N30_f0.txt");
time0 = table2array(time0);
time0 = 0.005.*(time0);
t0 = mean(time0);
time0_1 = readtable("N30_f0_1.txt");
time0_1 = table2array(time0_1);
time0_1 = 0.005.*(time0_1);
t0_1 = mean(time0_1);
time0_2 = readtable("N30_f0_2.txt");
time0_2 = table2array(time0_2);
time0_2 = 0.005.*(time0_2);
t0_2 = mean(time0_2);
time0_3 = readtable("N30_f0_3.txt");
time0_3 = table2array(time0_3);
time0_3 = 0.005.*(time0_3);
t0_3 = mean(time0_3);
time0_4 = readtable("N30_f0_4.txt");
time0_4 = table2array(time0_4);
time0_4 = 0.005.*(time0_4);
t0_4 = mean(time0_4);


times30 = [time0, time0_1, time0_2, time0_3, time0_4];
TT30 = [t0, t0_1, t0_2, t0_3, t0_4];
t0_on30 = t0;
TT30 = TT30.*5;
FF30 = [0, 0.1, 0.2, 0.3, 0.4];


% N = 40 (L/l_P = 8)
time0 = readtable("N40_f0.txt");
time0 = table2array(time0);
time0 = 0.005.*(time0);
t0 = mean(time0);
time0_1 = readtable("N40_f0_1.txt");
time0_1 = table2array(time0_1);
time0_1 = 0.005.*(time0_1);
t0_1 = mean(time0_1);
time0_2 = readtable("N40_f0_2.txt");
time0_2 = table2array(time0_2);
time0_2 = 0.005.*(time0_2);
t0_2 = mean(time0_2);
time0_3 = readtable("N40_f0_3.txt");
time0_3 = table2array(time0_3);
time0_3 = 0.005.*(time0_3);
t0_3 = mean(time0_3);


times40 = [time0, time0_1, time0_2, time0_3];
TT40 = [t0, t0_1, t0_2, t0_3]; 
t0_on40 = t0;
TT40 = TT40.*20;
FF40 = [0, 0.1, 0.2, 0.3];


%% %%%%%%%%% BLUMBERG (two-state) %%%%%%%%%%%%%%%
subplot(1,3,1)
hold on

% L/l_P = 2.2
L = 11;                 %polymer length      

ext = @(f) L.*((4/3) - 4 ./ (3 .* sqrt(f .* l_P .* beta + 1)) - ...
    (10 .* exp((900 ./ (beta .* f .* l_P)).^(1./4))) ./ ...
    (sqrt(f .* l_P .* beta) .* (exp((900 ./ (beta .* f .* l_P)).^(1./4)) - 1).^2) + ...
    ((f .* l_P .* beta).^1.62) ./ (3.55 + 3.8 .* (f .* l_P .* beta).^2.2));                     %force-extension

betaF = @(x) 1./l_P.*(-L.^2./(4.*(x - L)) - x./(4) + x.^2./(2.*L));                             %free energy F without force evaluated in x
betaG = @(f) betaF(ext(f)) - beta.* ext(f) .* f;                                                %free energy G with force evaluated in f
betaG1 = @(f) betaF(1) - beta.* 1.* f;                                                          %free energy G evaluated at f(x = 1)
DeltabetaG = @(f) -betaG(f) + betaG1(f);                                                        %free energy difference

tau = @(f) exp(DeltabetaG(f));                                                                  %looping time
tau = @(f) tau(f)./tau(FF11(end)).*TT11(end);                                                   %looping time rescaled to data point at highest force
fplot(tau, 'LineWidth', 2, 'Color', lightblue);

hold on

% L/l_P = 4.6
L = 23; 

ext = @(f) L.*((4/3) - 4 ./ (3 .* sqrt(f .* l_P .* beta + 1)) - ...
    (10 .* exp((900 ./ (beta .* f .* l_P)).^(1./4))) ./ ...
    (sqrt(f .* l_P .* beta) .* (exp((900 ./ (beta .* f .* l_P)).^(1./4)) - 1).^2) + ...
    ((f .* l_P .* beta).^1.62) ./ (3.55 + 3.8 .* (f .* l_P .* beta).^2.2));

betaF = @(x) 1./l_P.*(-L.^2./(4.*(x - L)) - x./(4) + x.^2./(2.*L));
betaG = @(f) betaF(ext(f)) - beta.* ext(f) .* f;
betaG1 = @(f) betaF(1) - beta.* 1.* f;
DeltabetaG = @(f) -betaG(f) + betaG1(f);

tau = @(f) exp(DeltabetaG(f));
tau = @(f) tau(f)./tau(FF23(end)).*TT23(end);
fplot(tau, 'LineWidth', 2, 'Color', lightblue);

% L/l_P = 6
L = 30;

ext = @(f) L.*((4/3) - 4 ./ (3 .* sqrt(f .* l_P .* beta + 1)) - ...
    (10 .* exp((900 ./ (beta .* f .* l_P)).^(1./4))) ./ ...
    (sqrt(f .* l_P .* beta) .* (exp((900 ./ (beta .* f .* l_P)).^(1./4)) - 1).^2) + ...
    ((f .* l_P .* beta).^1.62) ./ (3.55 + 3.8 .* (f .* l_P .* beta).^2.2));

betaF = @(x) 1./l_P.*(-L.^2./(4.*(x - L)) - x./(4) + x.^2./(2.*L));
betaG = @(f) betaF(ext(f)) - beta.* ext(f) .* f;
betaG1 = @(f) betaF(1) - beta.* 1.* f;
DeltabetaG = @(f) -betaG(f) + betaG1(f);

tau = @(f) exp(DeltabetaG(f));
tau = @(f) tau(f)./tau(FF30(end)).*TT30(end);
fplot(tau, 'LineWidth', 2, 'Color', lightblue);

% L/l_P = 8
L = 40;

ext = @(f) L.*((4/3) - 4 ./ (3 .* sqrt(f .* l_P .* beta + 1)) - ...
    (10 .* exp((900 ./ (beta .* f .* l_P)).^(1./4))) ./ ...
    (sqrt(f .* l_P .* beta) .* (exp((900 ./ (beta .* f .* l_P)).^(1./4)) - 1).^2) + ...
    ((f .* l_P .* beta).^1.62) ./ (3.55 + 3.8 .* (f .* l_P .* beta).^2.2));

betaF = @(x) 1./l_P.*(-L.^2./(4.*(x - L)) - x./(4) + x.^2./(2.*L));
betaG = @(f) betaF(ext(f)) - beta.* ext(f) .* f;
betaG1 = @(f) betaF(1) - beta.* 1.* f;
DeltabetaG = @(f) -betaG(f) + betaG1(f);

tau = @(f) exp(DeltabetaG(f));
tau = @(f) tau(f)./tau(FF40(end)).*TT40(end);
fplot(tau, 'LineWidth', 2, 'Color', lightblue);

% define errors based on standard error
err11 = std(times11)./sqrt(100);
err23 = std(times23)./sqrt(100).*2;
err30 = std(times30)./sqrt(100).*5;
err40 = std(times40)./sqrt(100).*20;

% plot error bars
errorbar(FF11,TT11,err11,"LineStyle","none", 'Marker','o','LineWidth',1,'LineStyle','none','Color',[0 0 0])
hold on
errorbar(FF23,TT23,err23,"LineStyle","none", 'Marker','^','LineWidth',1,'LineStyle','none','Color',[0 0 0])
errorbar(FF30,TT30,err30,"LineStyle","none", 'Marker','s','LineWidth',1,'LineStyle','none','Color',[0 0 0])
errorbar(FF40,TT40,err40,"LineStyle","none", 'Marker','d','LineWidth',1,'LineStyle','none','Color',[0 0 0])

xlim([0 1])
ylim([3000 2*10^7])
set(gca, 'YScale', 'log');

lgnd = legend({'', '', '', '', '$L/l_\mathrm{P} = 2.2$', '$L/l_\mathrm{P} = 4.6$', '$L/l_\mathrm{P} = 6$', '$L/l_\mathrm{P} = 8$'}, ...
              'Location','southeast','Interpreter','latex','FontSize',14);
set(lgnd,'color','none');

box on
ax = gca;
ax.FontSize = 14;

xlabel('$f$ [fN]', 'Interpreter','latex', 'FontSize',18)
ylabel('$\tau(f) \times C$ [s]', 'Interpreter','latex', 'FontSize',18)
title('Two-state', 'Interpreter','latex', 'FontSize',24)

% rescale force values (force of 1 in dimensionless units corresponds to 400 fN
xticks(0:0.2:1)
xt = get(gca, 'XTick');
set(gca,'XTickLabel',xt*400)

% rescale time values to match the data of Physical Review Letters 104, 258103 (2010)
% this predicts for L/l_P = 2.2 at f = 80 fN a looping time of 5.1 s
ytick_vals = [10^1 10^2 10^3 10^4 10^5] * TT11(3)/5.1;
yticks(ytick_vals);
ytick_scaled = ytick_vals / TT11(3) * 5.1;
ytick_labels = arrayfun(@(x) sprintf('$10^{%d}$', round(log10(x))), ytick_scaled, 'UniformOutput', false);

set(gca, 'YTickLabel', ytick_labels);
set(gca, 'TickLabelInterpreter', 'latex'); 



%% %%%%%%%%% SHIN (barrier-escape) %%%%%%%%%%%%%%%
subplot(1,3,2)

% L/l_P = 2.2
L = 11;
f_range = [0:0.01:1.2];
Lcut = 0;
x_c = 1;

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

    F = @(r) -log(r.^2.*(1 - (r./L).^2).^(-9./2).*exp(-3.*L./(4.*l_P.*(1 - (r./L).^2)))) - beta.*f.*r;


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

end

hold on
plot(f_range, tau./tau(101).*TT11(end), 'LineWidth', 2, 'Color', lightbrown)

% L/l_P = 4.6
L = 23;
f_range = [0:0.01:0.8];

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

    F = @(r) -log(r.^2.*(1 - (r./L).^2).^(-9./2).*exp(-3.*L./(4.*l_P.*(1 - (r./L).^2)))) - beta.*f.*r;


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

end

hold on
plot(f_range, tau./tau(61).*TT23(end), 'LineWidth', 2, 'Color', lightbrown)


% L/l_P = 6;
L = 30;
f_range = [0:0.01:0.6];

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

    F = @(r) -log(r.^2.*(1 - (r./L).^2).^(-9./2).*exp(-3.*L./(4.*l_P.*(1 - (r./L).^2)))) - beta.*f.*r;


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

end

hold on
plot(f_range, tau./tau(41).*TT30(end), 'LineWidth', 2, 'Color', lightbrown)


% L/l_P = 8
L = 40;
f_range = [0:0.01:0.5];

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

    F = @(r) -log(r.^2.*(1 - (r./L).^2).^(-9./2).*exp(-3.*L./(4.*l_P.*(1 - (r./L).^2)))) - beta.*f.*r;


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

end

hold on
plot(f_range, tau./tau(31).*TT40(end), 'LineWidth', 2, 'Color', lightbrown)

% plot erros
errorbar(FF11,TT11,err11,"LineStyle","none", 'Marker','o','LineWidth',1,'LineStyle','none','Color',[0 0 0])
hold on
errorbar(FF23,TT23,err23,"LineStyle","none", 'Marker','^','LineWidth',1,'LineStyle','none','Color',[0 0 0])
errorbar(FF30,TT30,err30,"LineStyle","none", 'Marker','s','LineWidth',1,'LineStyle','none','Color',[0 0 0])
errorbar(FF40,TT40,err40,"LineStyle","none", 'Marker','d','LineWidth',1,'LineStyle','none','Color',[0 0 0])

xlim([0 1])
ylim([3000 2*10^7])
set(gca, 'YScale', 'log');

lgnd = legend({'', '', '', '', '$L/l_\mathrm{P} = 2.2$', '$L/l_\mathrm{P} = 4.6$', '$L/l_\mathrm{P} = 6$', '$L/l_\mathrm{P} = 8$'}, ...
              'Location','southeast','Interpreter','latex','FontSize',14);

box on
ax = gca;
ax.FontSize = 14; 

xlabel('$f$ [fN]', 'Interpreter','latex', 'FontSize',18)
ylabel('$\tau(f) \times C$ [s]', 'Interpreter','latex', 'FontSize',18)

set(lgnd,'color','none');

title('Barrier escape', 'Interpreter','latex', 'FontSize',24)

% rescale force values to real units
xticks(0:0.2:1)
xt = get(gca, 'XTick');
set(gca,'XTickLabel',xt*400)


% rescale time values to real units
ytick_vals = [10^1 10^2 10^3 10^4 10^5] * TT11(3)/5.1;
yticks(ytick_vals);
ytick_scaled = ytick_vals / TT11(3) * 5.1;
ytick_labels = arrayfun(@(x) sprintf('$10^{%d}$', round(log10(x))), ytick_scaled, 'UniformOutput', false);

set(gca, 'YTickLabel', ytick_labels);
set(gca, 'TickLabelInterpreter', 'latex');  % This keeps font consistent across all labels


%% %%%%%%%%% LAEREMANS (inverse probability scaling) %%%%%%%%%%%%%%%
subplot(1,3,3)
hold on

% L/l_P = 2.2
L = 11;
ext = @(f) L.*((4/3) - 4 ./ (3 .* sqrt(f .* l_P .* beta + 1)) - ...
    (10 .* exp((900 ./ (beta .* f .* l_P)).^(1./4))) ./ ...
    (sqrt(f .* l_P .* beta) .* (exp((900 ./ (beta .* f .* l_P)).^(1./4)) - 1).^2) + ...
    ((f .* l_P .* beta).^1.62) ./ (3.55 + 3.8 .* (f .* l_P .* beta).^2.2));
betaF = @(x) 1./l_P.*(-L.^2./(4.*(x - L)) - x./(4) + x.^2./(2.*L));
betaG = @(f) betaF(ext(f)) - beta.* ext(f) .* f;
DeltabetaG = @(f) -betaG(f);
tau = @(f) exp(DeltabetaG(f));
tau = @(f) tau(f)./tau(FF11(end)).*TT11(end);
fplot(tau, 'LineWidth', 2, 'Color', purple);

hold on

% L/l_P = 3.6
L = 23;
ext = @(f) L.*((4/3) - 4 ./ (3 .* sqrt(f .* l_P .* beta + 1)) - ...
    (10 .* exp((900 ./ (beta .* f .* l_P)).^(1./4))) ./ ...
    (sqrt(f .* l_P .* beta) .* (exp((900 ./ (beta .* f .* l_P)).^(1./4)) - 1).^2) + ...
    ((f .* l_P .* beta).^1.62) ./ (3.55 + 3.8 .* (f .* l_P .* beta).^2.2));
betaF = @(x) 1./l_P.*(-L.^2./(4.*(x - L)) - x./(4) + x.^2./(2.*L));
betaG = @(f) betaF(ext(f)) - beta.* ext(f) .* f;
DeltabetaG = @(f) -betaG(f);
tau = @(f) exp(DeltabetaG(f));
tau = @(f) tau(f)./tau(FF23(end)).*TT23(end);
fplot(tau, 'LineWidth', 2, 'Color', purple);

% L/l_P = 6
L = 30;
ext = @(f) L.*((4/3) - 4 ./ (3 .* sqrt(f .* l_P .* beta + 1)) - ...
    (10 .* exp((900 ./ (beta .* f .* l_P)).^(1./4))) ./ ...
    (sqrt(f .* l_P .* beta) .* (exp((900 ./ (beta .* f .* l_P)).^(1./4)) - 1).^2) + ...
    ((f .* l_P .* beta).^1.62) ./ (3.55 + 3.8 .* (f .* l_P .* beta).^2.2));
betaF = @(x) 1./l_P.*(-L.^2./(4.*(x - L)) - x./(4) + x.^2./(2.*L));
betaG = @(f) betaF(ext(f)) - beta.* ext(f) .* f;
DeltabetaG = @(f) -betaG(f);
tau = @(f) exp(DeltabetaG(f));
tau = @(f) tau(f)./tau(FF30(end)).*TT30(end);
fplot(tau, 'LineWidth', 2, 'Color', purple);

% L/l_P = 8
L = 40;
ext = @(f) L.*((4/3) - 4 ./ (3 .* sqrt(f .* l_P .* beta + 1)) - ...
    (10 .* exp((900 ./ (beta .* f .* l_P)).^(1./4))) ./ ...
    (sqrt(f .* l_P .* beta) .* (exp((900 ./ (beta .* f .* l_P)).^(1./4)) - 1).^2) + ...
    ((f .* l_P .* beta).^1.62) ./ (3.55 + 3.8 .* (f .* l_P .* beta).^2.2));
betaF = @(x) 1./l_P.*(-L.^2./(4.*(x - L)) - x./(4) + x.^2./(2.*L));
betaG = @(f) betaF(ext(f)) - beta.* ext(f) .* f;
DeltabetaG = @(f) -betaG(f);
tau = @(f) exp(DeltabetaG(f));
tau = @(f) tau(f)./tau(FF40(end)).*TT40(end);
fplot(tau, 'LineWidth', 2, 'Color', purple);

% plot errors
errorbar(FF11,TT11,err11,"LineStyle","none", 'Marker','o','LineWidth',1,'LineStyle','none','Color',[0 0 0])
hold on
errorbar(FF23,TT23,err23,"LineStyle","none", 'Marker','^','LineWidth',1,'LineStyle','none','Color',[0 0 0])
errorbar(FF30,TT30,err30,"LineStyle","none", 'Marker','s','LineWidth',1,'LineStyle','none','Color',[0 0 0])
errorbar(FF40,TT40,err40,"LineStyle","none", 'Marker','d','LineWidth',1,'LineStyle','none','Color',[0 0 0])

xlim([0 1])
ylim([3000 2*10^7])
set(gca, 'YScale', 'log');


lgnd = legend({'', '', '', '', '$L/l_\mathrm{P} = 2.2$', '$L/l_\mathrm{P} = 4.6$', '$L/l_\mathrm{P} = 6$', '$L/l_\mathrm{P} = 8$'}, ...
              'Location','southeast','Interpreter','latex','FontSize',14);

box on
ax = gca;
ax.FontSize = 14; 

xlabel('$f$ [fN]', 'Interpreter','latex', 'FontSize',18)
ylabel('$\tau(f) \times C$ [s]', 'Interpreter','latex', 'FontSize',18)

set(lgnd,'color','none');

title('Inverse probability', 'Interpreter','latex', 'FontSize',24)

% convert forces to real units
xticks(0:0.2:1)
xt = get(gca, 'XTick');
set(gca,'XTickLabel',xt*400)

% convert timescales to real units
ytick_vals = [10^1 10^2 10^3 10^4 10^5] * TT11(3)/5.1;
yticks(ytick_vals);
ytick_scaled = ytick_vals / TT11(3) * 5.1;
ytick_labels = arrayfun(@(x) sprintf('$10^{%d}$', round(log10(x))), ytick_scaled, 'UniformOutput', false);

set(gca, 'YTickLabel', ytick_labels);
set(gca, 'TickLabelInterpreter', 'latex');  % This keeps font consistent across all labels

set(hfig, 'Position', [1.8,205,1534.4,456.8])
