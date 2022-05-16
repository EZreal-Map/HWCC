clc;clear
close all
load P_peak.mat  % P-->概率 peak-->峰量
% load new_ab
% P = flipud(b);
% peak = flipud(a);
%% 粒子群算法PSO
%% 初始化
N = 100;     % 种群中粒子个数f
D = 2;       % 粒子维数
T = 200;     % 最大迭代次数
c1 = 1.5;    % 学习因子1
c2 = 1.5;    % 学习因子2
Wmax = 2;    % 惯性权重
Wmin = 0.4;
Xmax = [4 10];   % 位置最大值
Xmin = [0 0];    % 位置最小值
Vmax = 0.1*(Xmax-Xmin);   % 速度最大值
Vmin = -Vmax;             % 速度最小值

%% 初始化种群个体（位置和速度）
x = rand(N,D).*(Xmax-Xmin)+Xmin;
v = rand(N,D).*(Vmax-Vmin)+Vmin;

% 记录种群中个体历史最优位置，因为所有个体目前只有初始位置x
p = x;
% 记录种群个体最优位置对应的 fit 值
p_best = zeros(N,1); % 开辟空间
for i = 1:N
    Cv = p(i,1);
    Cs = p(i,2);
    p_best(i) = fit_fun(peak,P,Cv,Cs);
end

g = zeros(1,D); % 开辟空间
g_best = inf;   % 开辟空间
for i = 1:N
    if p_best(i) < g_best
        g = p(i,:); % 记录更好的位置
        g_best = p_best(i); % 记录更好位置对应的fit值
    end
end

g_T = zeros(T,D);     % 开辟空间，记录不同迭代次数T的g
gbest_T = zeros(T,1); % 开辟空间，记录不同迭代次数T的g_best

for i = 1:T
    w = Wmax-(Wmax-Wmin)*i/T;
    for j =1:N
        % 更新个体最优位置和最优值
        temp = fit_fun(peak,P,x(j,1),x(j,2));
        if temp < p_best(j)
            p(j,:) = x(j,:);
            p_best(j) = temp;
            % 更新全局最优位置和最优值
            if temp < g_best
                g = p(j,:);
                g_best = temp;
            end
        end
        % 更新位置和速度
        v(j,:) = w*v(j,:)+c1*rand*(p(j,:)-x(j,:))+c2*rand*(g-x(j,:));
        for k = 1:D
            if v(j,k) > Vmax(k)
                v(j,k) = Vmax(k);
            elseif v(j,k) < Vmin(k)
                v(j,k) = Vmin(k);
            end
        end
        
        x(j,:) = x(j,:)+v(j,:);
        for k = 1:D
            if x(j,k) > Xmax(k)
                x(j,k) = Xmax(k);
            elseif x(j,k) < Xmin(k)
                x(j,k) = Xmin(k);
            end
        end
    end
    g_T(i,:) = g;
    gbest_T(i) = g_best;
    
end
%% 结果展示
f1 = figure;
plot(gbest_T)
fprintf(" Cv=%.4f Cs=%.4f\n ",g(1),g(2))


f2 = figure;
plot(P,peak);
%% 用矩法公式计算初值
n = length(P);
x_bar = mean(peak);
K = peak/x_bar;
Cv0 = sqrt(sum((K-1).^2)/(n-1));
Cs0 = sum((K-1).^3)/((n-3)*Cv0^3);
[~,new_Y] = fit_fun(peak,P,Cv0,Cs0);
% phi_p = Cs0/2*gaminv(1-P,4/Cs0^2,1)-2/Cs0;
% new_peak = x_bar*(1+Cv0*phi_p);
hold on
plot(P,new_Y)
% obj_function(peak,new_peak,1)
% fit_fun(peak,P,Cv0,Cs0)
[~,new_Y] = fit_fun(peak,P,g(1),g(2));
plot(P,new_Y)
legend('实测流量频率曲线','距估计初始值适线','粒子群适线')
hold off
% f3 = figure;
% plot(P,peak)
% for i=1:10:T
%
%     hold on
%     [~,new_Y] = fit_fun(peak,P,g_T(i,1),g_T(i,2));
%     plot(P,new_Y)
%     title(sprintf("第%d次迭代",i))
%     pause(0.1)
%
% end

ax = plot_Hpg(P,peak);
hold on
plot_Hpg(ax,P,new_Y);

%% 函数区
function ERR = obj_function(y,new_y,switch_fun)
switch switch_fun
    case 1
        % ERR = |y-new_y| 一次 绝对误差
        ERR = sum(abs(y-new_y));
    case 2
        % ERR = |y-new_y|/y 一次 相对误差
        ERR = sum(abs(y-new_y)./y);
    case 3
        % ERR = (y-new_y)^2 二次 绝对误差
        ERR = sum((y-new_y).^2);
    case 4
        % ERR = (y-new_y)^2/y 二次 相对误差
        ERR = sum((y-new_y).^2./y);
end
end


function [fit,new_Y] = fit_fun(Y,P,Cv,Cs)
Y_bar = mean(Y);
phi_p = Cs/2*gaminv(1-P,4/Cs^2,1)-2/Cs;
new_Y = Y_bar*(1+Cv*phi_p);
fit = obj_function(Y,new_Y,1); % 用选定的目标函数求误差
end


function ax=plot_Hpg(ax,P,y)
line_style = '-';
if nargin == 2
    y=P;
    P=ax;
    figure;
    ax = axes;
    line_style = 'x';
end

new_P = norminv(P,0,1);
Up_01 = norminv(0.01/100,0,1);
new_P = new_P-Up_01;
plot(ax,new_P,y,line_style)

x = [0.01 0.05 0.1 0.5 1.0 5 10 15 20 30 40 50 60 70 80 85 90 95 99 99.5 99.9];
Up = norminv(x/100,0,1);
Up_01 = norminv(0.01/100,0,1);
Lp = Up-Up_01;
ax.YGrid = 'on';
ax.XLim = [0 Lp(end)];
ax.XTick = Lp;
ax.XTickLabel = x';
ax.XGrid = 'on';
ax.XMinorGrid = 'on';
end