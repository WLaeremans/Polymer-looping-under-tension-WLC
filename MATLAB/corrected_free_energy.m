%% %%%%%%%%%%%% FIG A ***********
clc;
close all;

% plot simulation data free energy f = 0, N = 23 * with fit
scatter(centers_N23_f0./23, F_r_N_23_f0, 50, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.8, 0.8, 1], 'LineWidth', 1.5);
hold on
h = plot(fit_N23_f0_poly9);
set(h, 'Color', [0.04 0.04 0.58]);  % RGB color

% plot simulation data free energy f = 0.3, N = 23 *
scatter(centers_N23_f0_3./23, F_r_N_23_f0_3, 50, '^', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1, 0.8, 0.8], 'LineWidth', 1.5);
hold on

% load fit for N = 23, f = 0 and find the force-dependence via Eq. 23 in the article
        %fit N = 23, f = 0
    p1 = fittedmodelN23_f0.p1;
    p2 = fittedmodelN23_f0.p2;
    p3 = fittedmodelN23_f0.p3;
    p4 = fittedmodelN23_f0.p4;
    p5 = fittedmodelN23_f0.p5;
    p6 = fittedmodelN23_f0.p6;
    p7 = fittedmodelN23_f0.p7;
    p8 = fittedmodelN23_f0.p8;
    p9 = fittedmodelN23_f0.p9;
    p10 = fittedmodelN23_f0.p10;

    fit0 = @(x) p1*x.^9 + p2*x.^8 + p3*x.^7 + p4*x.^6 + p5*x.^5 + p6*x.^4 + p7*x.^3 + p8*x.^2 + p9*x + p10;

F0_3 = @(r) fit0(r) - log(sinh(0.3.*r)./(0.3.*r));

% normalize the probability distribution P, where the free energy is beta*F = -ln(P)
r_val = 0.01:0.01:25;           % Discretize domain
P0_3 = exp(-F0_3(r_val'));
area = trapz(r_val, P0_3);      % Compute numerical integral using trapezoidal rule
P_norm = P0_3 ./ area;          % Normalize so total area is 1
F_norm = -log(P_norm);
h = plot(r_val./23, F_norm);
set(h, 'Color', [0.58 0.04 0.04]);  % RGB color

% plot f = 0.6, N = 23
scatter(centers_N23_f0_6./23, F_r_N_23_f0_6, 50, 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.8, 1, 0.8], 'LineWidth', 1.5);
hold on
F0_6 = @(r) fit0(r) - log(sinh(0.6.*r)./(0.6.*r));
r_val = 0.01:0.01:25;           % Discretize domain
P0_6 = exp(-F0_6(r_val'));
area = trapz(r_val, P0_6);      % Compute numerical integral using trapezoidal rule
P_norm = P0_6 ./ area;          % Normalize so total area is 1
F_norm = -log(P_norm);
h = plot(r_val./23, F_norm);
set(h, 'Color', [0.04 0.58 0.04]);  % RGB color


xline(1)                    %capture distance
xline(1/23)                 %maximal distance (normalized)

box on
ax = gca;
ax.FontSize = 14; 

xlabel('$r/L$', 'Interpreter','latex', 'FontSize',18)
ylabel('$\beta F(r)$', 'Interpreter','latex', 'FontSize',18)

xlim([0 1.1])
ylim([0 17])

lgnd = legend({'$f = 0$~fN', '', '$f = 120$~fN', '', '$f = 240$~fN', ''}, 'Position',[0.470998641923463 0.706917988010185 0.27764557329593 0.171904765719459],'Interpreter','latex','FontSize',14)
set(gca, 'TickLabelInterpreter', 'latex');  % This keeps font consistent across all labels

annotation('textbox',...
    [0.172428571428571 0.821904760269896 0.108333611601882 0.0923809540158227],...
    'String',{'(a)'},...
    'Interpreter','latex',...
    'FontSize',18,...
    'EdgeColor','none');

%% %%%%%%%%%%%% PART B *********


figure
hold on
% Save the current folder
oldFolder = pwd;

% Change to the 'sim_data_11' folder

cd('sim_data_11');

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


% Go back to the original folder
cd(oldFolder);
times11 = [time0, time0_1, time0_2, time0_3, time0_4, time0_5, time0_6, time0_7, time0_8, time0_9, time1];
TT11 = [t0, t0_1, t0_2, t0_3, t0_4, t0_5, t0_6, t0_7, t0_8, t0_9, t1]; 
t0_on11 = t0;
TT11 = TT11.*1;
FF11 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
err11 = std(times11)./sqrt(100);

% plot inverse scaling with looping probability
purple = [160/256, 156/256, 252/256];

L = 11;
l_P = 5;
beta = 1;

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


%% numerical integration for barrier escape with correct free energy N = 11

        %fit f = 0
    p1 = fittedmodelN11_f0.p1;
    p2 = fittedmodelN11_f0.p2;
    p3 = fittedmodelN11_f0.p3;
    p4 = fittedmodelN11_f0.p4;
    p5 = fittedmodelN11_f0.p5;
    p6 = fittedmodelN11_f0.p6;
    p7 = fittedmodelN11_f0.p7;
    p8 = fittedmodelN11_f0.p8;
    p9 = fittedmodelN11_f0.p9;
    p10 = fittedmodelN11_f0.p10;

fit0N11 = @(x) p1*x.^9 + p2*x.^8 + p3*x.^7 + p4*x.^6 + p5*x.^5 + p6*x.^4 + p7*x.^3 + p8*x.^2 + p9*x + p10;

f_range = [0.001, 0.02:0.02:1];
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

L = 11;
% Initialize the result
result = 0;

tau = zeros(length(f_range),1);

for f_index = 1:length(f_range)
    result = 0;

    f = f_range(f_index);

        F = @(r) fit0N11(r) - log(sinh(f.*r)./(f.*r));      %free energy of Eq. 23


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
lightbrown = [1, 0.59, 0.59];

% Plot the line only (without markers)
p = plot(f_range, tau./tau(end).*TT11(end), '--', ...
    'LineWidth', 2, ...
    'Color', lightbrown);

hold on
errorbar(FF11,TT11,err11,"LineStyle","none", 'Marker','o','LineWidth',1.5,'LineStyle','none','Color',[0 0 0])
hold on

%% N = 23

% Save the current folder
oldFolder = pwd;

% Change to the 'sim_data_23' folder

cd('sim_data_23');


% N = 23
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

% Go back to the original folder
cd(oldFolder);

times23 = [time0, time0_1, time0_2, time0_3, time0_4, time0_5, time0_6];
TT23 = [t0, t0_1, t0_2, t0_3, t0_4, t0_5, t0_6]; 
t0_on23 = t0;
TT23 = TT23.*2;
FF23 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6];

err23 = std(times23)./sqrt(100).*2;


% plot inverse scaling with looping probability
L = 23;
Lcut = 2.3;
l_P = 5;
beta = 1;
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

xticks(0:0.2:1)
xt = get(gca, 'XTick');
set(gca,'XTickLabel',xt*400)

% ---------------------------
% Custom YTicks with LaTeX-style 10^x labels:
ytick_vals = [10^1 10^2 10^3 10^4 10^5] * TT11(3)/5.1;
yticks(ytick_vals);

% Now format the labels as powers of 10 (in LaTeX), but match font
ytick_scaled = ytick_vals / TT11(3) * 5.1;
ytick_labels = arrayfun(@(x) sprintf('$10^{%d}$', round(log10(x))), ytick_scaled, 'UniformOutput', false);

set(gca, 'YTickLabel', ytick_labels);
set(gca, 'TickLabelInterpreter', 'latex');  % This keeps font consistent across all labels

%% numerical integration barrier escape N = 23
f_range = [0.001, 0.02:0.02:0.6];

        %fit f = 0
    p1 = fittedmodelN23_f0.p1;
    p2 = fittedmodelN23_f0.p2;
    p3 = fittedmodelN23_f0.p3;
    p4 = fittedmodelN23_f0.p4;
    p5 = fittedmodelN23_f0.p5;
    p6 = fittedmodelN23_f0.p6;
    p7 = fittedmodelN23_f0.p7;
    p8 = fittedmodelN23_f0.p8;
    p9 = fittedmodelN23_f0.p9;
    p10 = fittedmodelN23_f0.p10;

    fit0 = @(x) p1*x.^9 + p2*x.^8 + p3*x.^7 + p4*x.^6 + p5*x.^5 + p6*x.^4 + p7*x.^3 + p8*x.^2 + p9*x + p10;


% Define the limits of integration
x_lower = x_c;
x_upper = L+Lcut;

% Define the function for the upper limit of y in terms of x
y_upper_func = @(x) x;

% Define the function for the upper limit of z in terms of y
z_upper_func = @(y) L+Lcut;

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

    F = @(r) fit0(r) - log(sinh(f.*r)./(f.*r));


    % Define the function to integrate over x, y, and z
    fun = @(x, y, z) exp(-(F(x) - F(y) + F(z)));

    integrand_Z = @(x) exp(-F(x));
    Z = integral(integrand_Z, x_c, L+Lcut);


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

lightbrown = [1, 0.59, 0.59];
% Plot the line only (without markers)
p = plot(f_range, tau./tau(end).*TT23(end), '--', ...
    'LineWidth', 2, ...
    'Color', lightbrown);
hold on

hold on
errorbar(FF23,TT23,err23,"LineStyle","none", 'Marker','^','LineWidth',1.5,'LineStyle','none','Color',[0 0 0])

xlim([0 1])
ylim([1000 2*10^7])
set(gca, 'YScale', 'log');

box on
ax = gca;
ax.FontSize = 14; 

xlabel('$f$ [fN]', 'Interpreter','latex', 'FontSize',18)
ylabel('$\tau(f) \times C$ [s]', 'Interpreter','latex', 'FontSize',18)

xticks(0:0.2:1)
xt = get(gca, 'XTick');
set(gca,'XTickLabel',xt*400)
xlim([0 1])

% ---------------------------
% Custom YTicks with LaTeX-style 10^x labels:
ytick_vals = [10^1 10^2 10^3 10^4 10^5] * TT11(3)/5.1;
yticks(ytick_vals);

% Now format the labels as powers of 10 (in LaTeX), but match font
ytick_scaled = ytick_vals / TT11(3) * 5.1;
ytick_labels = arrayfun(@(x) sprintf('$10^{%d}$', round(log10(x))), ytick_scaled, 'UniformOutput', false);

set(gca, 'YTickLabel', ytick_labels);
set(gca, 'TickLabelInterpreter', 'latex');  % This keeps font consistent across all labels


lgnd = legend({'Inv. prob.', 'Barr. escp.', 'Sim. $L/l_\mathrm{P} = 2.2$', '', '', 'Sim. $L/l_\mathrm{P} = 4.6$'}, ...
              'Location','southeast','Interpreter','latex','FontSize',14);
set(lgnd,'color','none');

annotation('textbox',...
    [0.172428571428571 0.821904760269896 0.108333611601882 0.0923809540158227],...
    'String',{'(b)'},...
    'Interpreter','latex',...
    'FontSize',18,...
    'EdgeColor','none');
