clc;clear
y_max = 10000;
y_gap = 100;
x = [0.01 0.05 0.1 0.5 1.0 5 10 15 20 30 40 50 60 70 80 85 90 95 99 99.5 99.9];
Up = norminv(x/100,0,1);
Up_01 = norminv(0.01/100,0,1);
Lp = Up-Up_01;

ax = axes;
ax.YGrid = 'on';
ax.XLim = [0 Lp(end)];
ax.XTick = Lp;
ax.XTickLabel = x';
ax.XGrid = 'on';
ax.XMinorGrid = 'on'

