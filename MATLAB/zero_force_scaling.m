%% Inverse scaling without force
clear;
clc;
close all;

purple = [160/256, 156/256, 252/256];
dark_grey = [0.5, 0.5, 0.5];

figure(3)

% read data
time5 = readtable("L5.txt");
time5 = table2array(time5);
time5 = 0.005.*(time5);
t5 = mean(time5);

time10 = readtable("L10.txt");
time10 = table2array(time10);
time10 = 0.005.*(time10);
time10(101:end) = [];
t10 = mean(time10);

time15 = readtable("L15.txt");
time15 = table2array(time15);
time15 = 0.005.*(time15);
t15 = mean(time15);

time20 = readtable("L20.txt");
time20 = table2array(time20);
time20 = 0.005.*(time20);
t20 = mean(time20);

time30 = readtable("L30.txt");
time30 = table2array(time30);
time30 = 0.005.*(time30);
t30 = mean(time30);

time50 = readtable("L50.txt");
time50 = table2array(time50);
time50 = 0.005.*(time50);
t50 = mean(time50);

time100 = readtable("L100.txt");
time100 = table2array(time100);
time100 = 0.005.*(time100);
t100 = mean(time100);

timesf0all = [time5, time10, time15, time20, time30, time50, time100];

errf0 = std(timesf0all)./sqrt(100);


Ld = [5, 10, 15, 20, 30, 50, 100];
tt = [t5, t10, t15,t20, t30, t50, t100];


hold on;
l_P = 5;

% evaluate inverse probability via algoritm of 
% Physics Letters A 381, 1029 (2017)
lmax = 10;
h = 0.005;
final = [];

for k = 0:0.1:80
    beta = 0.25 * k + 1;
    if beta <= 3
        Nmax = 90000;
    else
        Nmax = 9000;
    end
    L = [];

    for n = 0:Nmax
        f = h * n * 1i;
        H = zeros(lmax + 1);
        for i = 0:lmax
            for j = 0:lmax
                switch i - j
                    case -1
                        H(i + 1, j + 1) = f * (i + 1) / sqrt((2 * i + 1) * (2 * i + 3));
                    case 0
                        H(i + 1, j + 1) = i * (i + 1) / 2;
                    case 1
                        H(i + 1, j + 1) = f * i / sqrt((2 * i - 1) * (2 * i + 1));
                end
            end
        end
        M = expm(-beta * H);
        Z = M(1, 1);
        L = [L; Z];
    end
    L = real(L);
    Pz = [];
    P1z = [];

    for l = -2:1199
        xi = 0.001 * l;
        P = (h * beta / pi) * sum(L .* cos((0:Nmax)' * h * xi * beta));
        Pz = [Pz; xi, P];
        P1z = [P1z; P];
    end

    V = P1z;
    L1 = V(3:end);
    L2 = V(1:end-2);
    LPR = (L1 - L2) / (0.001 * 2);
    LPR = LPR(2:end);
    QR = zeros(length(LPR), 1);

    for i = 2:length(LPR) + 1
        QR(i - 1) = LPR(i - 1) * 1 / ((i - 1) * 0.001) * (-1 / (2 * pi));
    end

    QR1 = [(2:length(LPR) + 1)' * 0.001, LPR .* (1 ./ (2:length(LPR) + 1)') .* 1 ./ ((2:length(LPR) + 1)' * 0.001) * (-1 / (2 * pi))];

    final = [final; beta, (1 / beta^3) * QR(2)];
    k
end


errorbar(Ld,tt,errf0,"LineStyle","none", 'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
hold on
plot(final(11:end,1).*5, 12./final(11:end,2), 'LineWidth', 3, 'LineStyle','-', 'Color', purple);

% add approximated looping probability for small L of 
% Macromolecules 17, 689 (1984)
g = @(x) 0.19*(112.04*l_P^2/x^5*exp(-14.055*l_P/x + 0.246*x/l_P))^(-1); 
fplot(g, 'LineWidth', 2, 'LineStyle', '--', 'Color', dark_grey)
errorbar(Ld,tt,errf0,"LineStyle","none", 'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])


box on
ax = gca;
ax.FontSize = 14; 

xlabel('$L = Nb$ [LJ]', 'Interpreter','latex', 'FontSize',18)
ylabel('$\tau(f = 0)$ [LJ]', 'Interpreter','latex', 'FontSize',18)

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlim([0 100])
ylim([2*10^3 10^5])

% FJC result valid for long L
line([30,90],[30*30^(3/2),30*90^(3/2)], "LineWidth", 2, 'Color', 'k')
annotation('textbox',...
    [0.665642857142853 0.425014153021068 0.140171947670677 0.0792715612646545],...
    'String',{'$\sim L^{3/2}$'},...
    'Rotation',40,...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');

lgnd = legend({'Simulation', 'Inverse probability num.', 'Inverse probability an.'},'Location', 'northeast','Interpreter','latex','FontSize',14)
set(lgnd,'color','none');
