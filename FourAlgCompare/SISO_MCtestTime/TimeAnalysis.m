%*************************************************************************%
%   分数阶粒子滤波仿真复现
%   论文：     fractional order PF
%   目的：分数阶粒子滤波算法测试
%
%   函数实验:    D^{0.7} x_k = 3*sin(2*x_{k-1}) - exp(x_{k-1}) + w_k
%                        y_k = x_k + v_k
%   结果：较好的对状态进行估计
%
%   备注： 计算四种算法的运行时间及误差范数
%           
%*************************************************************************%

clc
clear

