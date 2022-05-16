clc;clear
x=[0.01:0.01:.1,0.2:.1:1,1.2:.2:2,3:20,22:2:80,81:98,98.2:.2:99,99.1:.1:99.9]; %横坐标的值
y=norminv(x/100,0,1);
y=y-y(1);
for i=1:size(y,2)
    line([y(i),y(i)],[0 19]);
    pause(0.01)
end

for i=1:19 % 10为纵坐标的格数，可以自已设置。
    line([0 y(end)],[i i]);
    pause(0.01)
end

h=findobj('type','axes');set(h,'xtick',[],'ytick',[],'xlim',[0 y(end)],'ylim',[0 ,19]);
xx=[0.01 0.05 0.1 0.5 1.0 5 10 15 20 30 40 50 60 70 80 85 90 95 99 99.5 99.9];%标横坐标的值，可以自己设置
yy=norminv(xx/100,0,1);yy=yy-yy(1);
for i=1:size(xx,2)
    text('string',num2str(xx(i)),'HorizontalAlignment','center','pos',[yy(i),-0.5]);
end