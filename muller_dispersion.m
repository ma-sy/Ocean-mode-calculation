function [f_and_c] = muller_dispersion(cw,rou1,alphaw,cp0,cs0,rou2,alphap,alphas,h,f_range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%变量说明%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cw：水中声速，单位m/s
% rou1：水密度，单位kg/m^3
% alphaw：水中声吸收系数，单位dB/波长
% cp0：海底压缩波速，单位m/s
% cs0：海底剪切波速，单位m/s
% rou2：海底密度，单位kg/m^3
% alphap：海底压缩波声吸收系数，单位dB/波长
% alphas：海底剪切波声吸收系数，单位dB/波长
% h：海深，单位m
% f_range：频率计算范围，单位Hz，以向量形式给出
% f_and_c：输出结果矩阵，频率和相速度，按行排列
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% 哈工程 马嗣宇 2022.9 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 计入声吸收衰减，引入复声速
    c1=cw/(1+1i*(alphaw*log(10)/40/pi));
    cs=cs0/(1+1i*(alphas*log(10)/40/pi));
    cp=cp0/(1+1i*(alphap*log(10)/40/pi));
    
    f_and_c=[];     % 输出结果矩阵，频率和相速度，按行排列
    delt=0.5;       % 局部峰值搜索复空间划分间隔
    ybxl=0.0001;    % 局部峰值搜索判定门限
    for f1=f_range
       %% 方程定义
        f0 = @(c) tanh(2.*pi.*f1.*h./c.*sqrt(1-c.^2/c1.^2)) .*rou1/rou2.*c.^4/cs0.^4 .*sqrt(1-c.^2/cp0.^2)...
            +(2-c.^2/cs0.^2).^2.*sqrt(1-c.^2/c1.^2)...
            -4.*sqrt(1-c.^2/cs0.^2).*sqrt(1-c.^2/cp0.^2).*sqrt(1-c.^2/c1.^2);
        g0 = @(c) 1./(1+abs(f0(c)));
        %% 判断范围
        c0test = [];
        %% 局部峰值搜索法
        for im=0:delt:10
            for re=0:delt:cs0
                if g0(1i*im+re)==max([g0(1i*(im-delt)+(re-delt)), g0(1i*(im-delt)+(re)), g0(1i*(im-delt)+(re+delt)), ...
                        g0(1i*(im)+(re-delt)),      g0(1i*im+re),          g0(1i*(im)+(re+delt)), ...
                        g0(1i*(im+delt)+(re-delt)), g0(1i*(im+delt)+(re)), g0(1i*(im+delt)+(re+delt))])
                    if (g0(1i*im+re)-(g0(1i*(im-delt)+(re-delt))+ g0(1i*(im-delt)+(re))+ g0(1i*(im-delt)+(re+delt))+ ...
                            g0(1i*(im)+(re-delt))+                            g0(1i*(im)+(re+delt))+ ...
                            g0(1i*(im+delt)+(re-delt))+ g0(1i*(im+delt)+(re))+ g0(1i*(im+delt)+(re+delt)))/8 )/g0(1i*im+re)>=ybxl
                        if re~=cw
                            c0test = [1i*im+re-delt 1i*im+re 1i*im+re+delt];
                           %% muller法
                            [res,~,~] = muller (f0,c0test);
                            f_and_c = [f_and_c,[f1;res]];
                        end
                    end
                end
            end
        end
    end
    figure;
    scatter(f_and_c(1,:),f_and_c(2,:),5,'filled');
    title({['h=',num2str(h),' c_w=',num2str(cw),' \rho_w=',num2str(rou1),' \alpha=',num2str(alphaw)];['c_p=',num2str(cp0),' c_s=',num2str(cs0),' \rho_e=',num2str(rou2),' \alpha_p=',num2str(alphap),' \alpha_s=',num2str(alphas)]});
    xlabel('f(Hz)');ylabel('c(m/s)')
end