clc;clear
close all
%% 数据的加载
filename='..\课程设计数据\No2019101629.xls';
cell_data=readcell(filename,'sheet','防洪计算');
% 水文站洪水特征统计 A3-D23
temp=cell_data(3:23,1:4);
colname_english={'year' 'peak' 'volume_1day' 'volume_3day'};
flood_characteristics=cell2struct(temp,colname_english,2); % 洪水特征struct
% flood_characteristics.size=size(flood_characteristics,1)
% 水文站典型洪水流量 F3-H32
temp=cell_data(3:32,6:8);
colname_english={'month_day' 'hour' 'flow'};
typical_flood=cell2struct(temp,colname_english,2); % 典型洪水struct

% 坝址有关数据 N1-N11
temp=cell_data(1:11,14);
rowname_english={'F' 'river_length' 'slope' 'longitude' ...
    'latitude' 'rainfall_point2area_1hour' 'rainfall_point2area_6hour' ...
    'control_area' 'M' 'blow' 'wind_speed'};
drainage_basin=cell2struct(temp,rowname_english,1); % 流域struct

% 特大洪水数据 M15-N16
temp=cell_data(15:16,13:14);
colname_english={'year' 'peak'};
catastrophic_flood=cell2struct(temp,colname_english,2); % 特大洪水struct

%% 数据的整理

%% 设计指导 一 1 （4）
temp_index=find([flood_characteristics.volume_3day]~=-1);
% corrcoef([flood_characteristics(temp_index).volume_1day],[flood_characteristics(temp_index).volume_3day])
% plot([flood_characteristics(temp_index).volume_1day],[flood_characteristics(temp_index).volume_3day],'.')
% 进行1日洪量与3日洪量的线性拟合（相关分析法）
[p,s]=polyfit([flood_characteristics(temp_index).volume_1day],[flood_characteristics(temp_index).volume_3day],1);

temp_index=find([flood_characteristics.volume_3day]==-1);
flood_characteristics(temp_index).volume_3day=polyval(p,[flood_characteristics(temp_index).volume_1day]);
% hold on
% plot([flood_characteristics(temp_index).volume_1day],[flood_characteristics(temp_index).volume_3day],'*r')
% hold off

%% 设计指导 一 2 
% 经验频率计算：下面计算步骤默认2个前提
% 1、历史大洪水只有一个调查期（不像书上有3个调查考证期）
% 2、实测洪水年份是连续的
%% 统一处理法 调查考证期 频率计算
a=size(catastrophic_flood,1); % N年中有a个特大值
peak=[catastrophic_flood.peak flood_characteristics.peak]';
year=[catastrophic_flood.year flood_characteristics.year]';
[peak,temp_index]=sort(peak,'descend');
year=year(temp_index);
N=max(year)-min(year)+1; % 最大年份与最小年份之差（注意不是最大洪水年份到最小洪水年份之差）

% 要求实测洪水年是连续的，一般都能满足，所以此处就默认按照一般情况处理（如不连续下面计算方法出错）
min_measured_year=min([flood_characteristics.year]);
max_measured_year=max([flood_characteristics.year]);
n=max_measured_year-min_measured_year+1; % 实测系列年数
l=0; %初始化l
for i=1:a
    if catastrophic_flood(i).year>=min_measured_year && catastrophic_flood(i).year <= max_measured_year
        l=l+1; % 来自于实测资料的特大洪水年个数
    end
end

Pm_1=zeros(n+a-l,1);
Pm_1(1:a)=(1:a)/(N+1); % 公式2-7
Pm_1(a+1:end)=Pm_1(a)+(1-Pm_1(a))*((l+1:n)-l)/(n-l+1); % 公式2-9

%% 分别处理法 调查考证期 频率计算
Pm_2=zeros(n+a-l,1);
Pm_2(1:a)=(1:a)/(N+1); % 公式2-7
Pm_2(a+1:end)=(l+1:n)/(n+1); % 公式2-10
% [Pm_1';Pm_2']
P=Pm_1';
peak=peak';
save('P_peak','P','peak')
