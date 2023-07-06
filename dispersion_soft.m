function [f_and_c] = dispersion_soft(c1,rou1,cp0,cs0,rou2,h,f_range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%变量说明%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c1：水中声速，单位m/s
% rou1：水密度，单位kg/m^3
% cp0：海底压缩波速，单位m/s
% cs0：海底剪切波速，单位m/s
% rou2：海底密度，单位kg/m^3
% h：海深，单位m
% f_range：频率计算范围，单位Hz，以向量形式给出
% f_and_c：输出结果矩阵，频率和相速度，按行排列
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% 哈工程 马嗣宇 2022.9 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    f_and_c=[];
    w_and_k=[];
    for f1=f_range
        %% 二分法
        ans_zs=[];
        a=0;
        b=cs0;
        f = @(c) tanh(2.*pi.*f1/c.*sqrt(1-c.^2/c1.^2).*h).*rou1/rou2.*c.^4/cs0.^4.*sqrt(1-c.^2/cp0.^2)+(2-c.^2/cs0.^2).^2.*sqrt(1-c.^2/c1.^2)-4.*sqrt(1-c.^2/cs0.^2).*sqrt(1-c.^2/cp0.^2).*sqrt(1-c.^2/c1.^2);
        for c=a:0.1:b
            ans_zs=[ans_zs,[c;f(c)]];
        end
        i=1;
        l=length(ans_zs(1,:));
        while i~=l-1
            if ans_zs(2,i).*ans_zs(2,i+1)<0
                a=ans_zs(1,i);
                b=ans_zs(1,i+1);
                t=(a+b)./2;
                while abs(b-a)>1e-8
                    if f(t)*f(b)<0
                        a=t;
                    else
                        b=t;
                    end
                    t=(a+b)./2;
                    c=t;
                end

                if abs(f(c)-0)<1e-5
                    f_and_c=[f_and_c,[f1;c]];
                    w_and_k=[w_and_k,[2*pi*f1;2*pi*f1/c]];
                end
            end
            i=i+1;
        end
    end
    figure;
%     scatter(f_and_c(1,:),f_and_c(2,:),5,'filled');
    plot(f_and_c(1,:),f_and_c(2,:),'b','LineWidth',1);
    hold on;
    plot(f_and_c(1,1:end-1),diff(w_and_k(1,:))./diff(w_and_k(2,:)),'r','LineWidth',1);
    title(['h=',num2str(h),' c_w=',num2str(c1),' \rho_w=',num2str(rou1),' c_p=',num2str(cp0),' c_s=',num2str(cs0),' \rho_e=',num2str(rou2),]);
    xlabel('f(Hz)');ylabel('c(m/s)')
end