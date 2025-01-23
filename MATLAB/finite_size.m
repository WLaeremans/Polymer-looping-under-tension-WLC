clear
clc
close all

% N = 11
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

times11 = [time0, time0_1, time0_2, time0_3, time0_4, time0_5, time0_6, time0_7, time0_8, time0_9, time1];
TT11 = [t0, t0_1, t0_2, t0_3, t0_4, t0_5, t0_6, t0_7, t0_8, t0_9, t1]; 
t0_on11 = t0;
TT11 = TT11./TT11(1);
FF11 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];

err11 = std(times11)./sqrt(100)./t0_on11;

errorbar(FF11,TT11,err11,"LineStyle","none", 'Marker','o','LineWidth',1,'LineStyle','none','Color',[0 0 0])

box on
ax = gca;
ax.FontSize = 12; 
xlabel('$\beta f b$', 'Interpreter','latex', 'FontSize',18)
ylabel('$\ln\left[ \tau_L(N,f)/\tau_L(N,f=0) \right]$', 'Interpreter','latex', 'FontSize',18)
set(gca, 'YScale', 'log');

% N = 22

hold on

time0fse = readtable("N22_f0_fse.txt");
time0fse = table2array(time0fse);
time0fse = 0.005.*(time0fse);
t0fse = mean(time0fse);

time0_1fse = readtable("N22_f01_fse.txt");
time0_1fse = table2array(time0_1fse);
time0_1fse = 0.005.*(time0_1fse);
t0_1fse = mean(time0_1fse);

time0_2fse = readtable("N22_f02_fse.txt");
time0_2fse = table2array(time0_2fse);
time0_2fse = 0.005.*(time0_2fse);
t0_2fse = mean(time0_2fse);

time0_3fse = readtable("N22_f03_fse.txt");
time0_3fse = table2array(time0_3fse);
time0_3fse = 0.005.*(time0_3fse);
t0_3fse = mean(time0_3fse);

time0_4fse = readtable("N22_f04_fse.txt");
time0_4fse = table2array(time0_4fse);
time0_4fse = 0.005.*(time0_4fse);
t0_4fse = mean(time0_4fse);

time0_5fse = readtable("N22_f05_fse.txt");
time0_5fse = table2array(time0_5fse);
time0_5fse = 0.005.*(time0_5fse);
t0_5fse = mean(time0_5fse);

time0_6fse = readtable("N22_f06_fse.txt");
time0_6fse = table2array(time0_6fse);
time0_6fse = 0.005.*(time0_6fse);
t0_6fse = mean(time0_6fse);

time0_7fse = readtable("N22_f07_fse.txt");
time0_7fse = table2array(time0_7fse);
time0_7fse = 0.005.*(time0_7fse);
t0_7fse = mean(time0_7fse);

time0_8fse = readtable("N22_f08_fse.txt");
time0_8fse = table2array(time0_8fse);
time0_8fse = 0.005.*(time0_8fse);
t0_8fse = mean(time0_8fse);

time0_9fse = readtable("N22_f09_fse.txt");
time0_9fse = table2array(time0_9fse);
time0_9fse = 0.005.*(time0_9fse);
t0_9fse = mean(time0_9fse);

time1fse = readtable("N22_f1_fse.txt");
time1fse = table2array(time1fse);
time1fse = 0.005.*(time1fse);
t1fse = mean(time1fse);

times22 = [time0fse, time0_1fse, time0_2fse, time0_3fse, time0_4fse, time0_5fse, time0_6fse, time0_7fse, time0_8fse, time0_9fse, time1fse];
TT22 = [t0, t0_1, t0_2, t0_3, t0_4, t0_5, t0_6, t0_7, t0_8, t0_9, t1]; 
t0_on22 = t0;
TT22 = TT22./TT22(1);
FF22 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];

err22 = std(times22)./sqrt(100)./t0_on22;

errorbar(FF22,TT22,err22,"LineStyle","none", 'Marker','^','LineWidth',1,'LineStyle','none','Color',[1 0 0])

set(gca, 'YScale', 'log');



lgnd = legend({'N = 11', 'N = 22'},'Location','southeast','Interpreter','latex','FontSize',14)

box on
ax = gca;
ax.FontSize = 14; 

xlabel('$f$ [LJ]', 'Interpreter','latex', 'FontSize',18)
ylabel('$\tau(f) / \tau(f=0)$', 'Interpreter','latex', 'FontSize',18)

set(lgnd,'color','none');

xlim([0 1])
ylim([0.7 400])
