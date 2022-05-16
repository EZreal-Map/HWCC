format = 'M月d日 HH:mm';
infmt = 'M.d H.mm'
datetime(2004,2,29,14,0,0,'Format',format);
t=datetime('5.12 23.21','InputFormat',infmt)
t.Minute=t.Minute*0.6
clc;clear
load('time.mat')
a = ans
a(1,:).Time