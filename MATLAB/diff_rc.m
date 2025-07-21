clear;
clc;
close all;

purple = [160/256, 156/256, 252/256];
darkbrown = [0.619607843137255 0.192156862745098 0.192156862745098];
lightbrown = [1, 0.59, 0.59];
lightblue = [0.62, 0.89, 1];

%% read data with rc = 1
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
err11 = std(times11)./sqrt(100)./t0_on11.*1;
hold on

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


times23 = [time0, time0_1, time0_2, time0_3, time0_4, time0_5, time0_6];
TT23 = [t0, t0_1, t0_2, t0_3, t0_4, t0_5, t0_6]; 
t0_on23 = t0;
TT23 = TT23./TT23(1).*2;
FF23 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
err23 = std(times23)./sqrt(100)./t0_on23.*2;

TT23rescaler = TT23;

%% plot
errorbar(FF11,TT11,err11,"LineStyle","none", 'Marker','o','LineWidth',1,'LineStyle','none','Color',[0 0 0])
errorbar(FF11, TT11, err11, 'LineWidth', 1, 'Color', [0 0 0])
errorbar(FF23,TT23,err23,"LineStyle","none", 'Marker','o','LineWidth',1,'LineStyle','none','Color',[0 0 0])
errorbar(FF23, TT23, err23, 'LineWidth', 1, 'Color', [0 0 0])

set(gca, 'YScale', 'log');

%% add data rc = 0.5
% N = 23
hold on
time0 = readtable("N23_f0_rc05.txt");
time0 = table2array(time0);
time0 = 0.005.*(time0);
t0rc = mean(time0);
time0_1 = readtable("N23_f0_1_rc05.txt");
time0_1 = table2array(time0_1);
time0_1 = 0.005.*(time0_1);
t0_1rc = mean(time0_1);
time0_2 = readtable("N23_f0_2_rc05.txt");
time0_2 = table2array(time0_2);
time0_2 = 0.005.*(time0_2);
t0_2rc = mean(time0_2);
time0_3 = readtable("N23_f0_3_rc05.txt");
time0_3 = table2array(time0_3);
time0_3 = 0.005.*(time0_3);
t0_3rc = mean(time0_3);
time0_4 = readtable("N23_f0_4_rc05.txt");
time0_4 = table2array(time0_4);
time0_4 = 0.005.*(time0_4);
t0_4rc = mean(time0_4);
time0_5 = readtable("N23_f0_5_rc05.txt");
time0_5 = table2array(time0_5);
time0_5 = 0.005.*(time0_5);
t0_5rc = mean(time0_5);
time0_6 = readtable("N23_f0_6_rc05.txt");
time0_6 = table2array(time0_6);
time0_6 = 0.005.*(time0_6);
t0_6rc = mean(time0_6);

times23 = [time0, time0_1, time0_2, time0_3, time0_4, time0_5, time0_6];
TT23 = [t0rc, t0_1rc, t0_2rc, t0_3rc, t0_4rc, t0_5rc, t0_6rc]; 
t0_on23 = t0rc;
TT23 = TT23./TT23(1).*2;
FF23 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
err23 = std(times23)./sqrt(100)./t0_on23.*2;
errorbar(FF23,TT23,err23,"LineStyle","none", 'Marker','^','LineWidth',1,'LineStyle','none','Color',[1 0 0])

errorbar(FF23, TT23, err23, 'LineWidth', 1, 'Color', [1 0 0])

% N = 11
hold on
time0 = readtable("N11_f0_rc05.txt");
time0 = table2array(time0);
time0 = 0.005.*(time0);
t0rc = mean(time0);
time0_1 = readtable("N11_f0_1_rc05.txt");
time0_1 = table2array(time0_1);
time0_1 = 0.005.*(time0_1);
t0_1rc = mean(time0_1);
time0_2 = readtable("N11_f0_2_rc05.txt");
time0_2 = table2array(time0_2);
time0_2 = 0.005.*(time0_2);
t0_2rc = mean(time0_2);
time0_3 = readtable("N11_f0_3_rc05.txt");
time0_3 = table2array(time0_3);
time0_3 = 0.005.*(time0_3);
t0_3rc = mean(time0_3);
time0_4 = readtable("N11_f0_4_rc05.txt");
time0_4 = table2array(time0_4);
time0_4 = 0.005.*(time0_4);
t0_4rc = mean(time0_4);
time0_5 = readtable("N11_f0_5_rc05.txt");
time0_5 = table2array(time0_5);
time0_5 = 0.005.*(time0_5);
t0_5rc = mean(time0_5);
time0_6 = readtable("N11_f0_6_rc05.txt");
time0_6 = table2array(time0_6);
time0_6 = 0.005.*(time0_6);
t0_6rc = mean(time0_6);
time0_7 = readtable("N11_f0_7_rc05.txt");
time0_7 = table2array(time0_7);
time0_7 = 0.005.*(time0_7);
t0_7rc = mean(time0_7);
time0_8 = readtable("N11_f0_8_rc05.txt");
time0_8 = table2array(time0_8);
time0_8 = 0.005.*(time0_8);
t0_8rc = mean(time0_8);
time0_9 = readtable("N11_f0_9_rc05.txt");
time0_9 = table2array(time0_9);
time0_9 = 0.005.*(time0_9);
t0_9rc = mean(time0_9);
time1 = readtable("N11_f1_rc05.txt");
time1 = table2array(time1);
time1 = 0.005.*(time1);
t1rc = mean(time1);

times11 = [time0, time0_1, time0_2, time0_3, time0_4, time0_5, time0_6, time0_7, time0_8, time0_9, time1];
TT11 = [t0rc, t0_1rc, t0_2rc, t0_3rc, t0_4rc, t0_5rc, t0_6rc, t0_7rc, t0_8rc, t0_9rc, t1rc]; 
t0_on11 = t0rc;
TT11 = TT11./TT11(1).*1;
FF11 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
err11 = std(times11)./sqrt(100)./t0_on11.*1;
errorbar(FF11,TT11,err11,"LineStyle","none", 'Marker','^','LineWidth',1,'LineStyle','none','Color',[1 0 0])
hold on
errorbar(FF11, TT11, err11, 'LineWidth', 1, 'Color', [1 0 0])

xlim([0 1])
ylim([0.7 4000])
set(gca, 'YScale', 'log');

lgnd = legend({'Sim. $r_\mathrm{c} = 1$', '','','','Sim. $r_\mathrm{c} = 0.5$'},'Location','southeast','Interpreter','latex','FontSize',14)

box on
ax = gca;
ax.FontSize = 14; 

xlabel('$f$ [LJ]', 'Interpreter','latex', 'FontSize',18)
ylabel('$\tau(f) / \tau(f=0) \times C$', 'Interpreter','latex', 'FontSize',18)

set(lgnd,'color','none');

annotation('textbox',...
    [0.305642857142857 0.627142855871293 0.201857489476893 0.0771428584144229],...
    'String',{'$L/l_\mathrm{P} = 4.6$'},...
    'Interpreter','latex',...
    'FontSize',14,...
    'EdgeColor','none');

annotation('textbox',...
    [0.657785714285713 0.407142855871295 0.201857489476893 0.0771428584144229],...
    'String',{'$L/l_\mathrm{P} = 2.2$'},...
    'Interpreter','latex',...
    'FontSize',14,...
    'EdgeColor','none');
