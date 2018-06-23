function [  ] = plotAperture( xValue, apertureValue)
% 画出口径面处的数值解
% 本函数用于供C++程序内部调用

%% 绘制口径面处的解
figure

plot(xValue,apertureValue,'b');


