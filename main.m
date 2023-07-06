clc;
close all;
clear all;

cw = input('请输入水中声速，单位m/s：');
rou1 = input('请输入水密度，单位kg/m3：');
alphaw = input('请输入水中声吸收系数，单位dB/波长：');
cp0 = input('请输入海底P波速，单位m/s：');
cs0 = input('请输入海底S波速，单位m/s：');
rou2 = input('请输海底密度，单位kg/m3：');
alphap = input('请输入海底P波声吸收系数，单位dB/波长：');
alphas = input('请输入海底S波声吸收系数，单位dB/波长：');
h = input('请输入海深，单位m：');
f_range = input('频率计算范围，单位Hz，如0:0.5:10：');
tic
if alphaw==0 && alphap==0 && alphas==0
    if cw<=cs0
        dispersion_hard(cw,rou1,cp0,cs0,rou2,h,f_range);
    else
        dispersion_soft(cw,rou1,cp0,cs0,rou2,h,f_range);
    end
else
    muller_dispersion(cw,rou1,alphaw,cp0,cs0,rou2,alphap,alphas,h,f_range);
end
toc