function [f_and_c] = dispersion_hard(c1,rou1,cp0,cs0,rou2,h,f_range)
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
    %% 二分法
    f_and_c=[];
    for f1=f_range
        %% 第一段
        ans_zs=[];
        a=0;
        b=c1;
        %% 方程定义，模态方程按硬质海底条件展开(第一段)
        f = @(c) tanh(2.*pi.*f1/c.*sqrt(1-c.^2/c1.^2).*h).*rou1/rou2.*c.^4/cs0.^4.*sqrt(1-c.^2/cp0.^2)+(2-c.^2/cs0.^2).^2.*sqrt(1-c.^2/c1.^2)-4.*sqrt(1-c.^2/cs0.^2).*sqrt(1-c.^2/cp0.^2).*sqrt(1-c.^2/c1.^2);
        for c=a:0.1:b
            ans_zs=[ans_zs,[c;f(c)]];
        end
        i=1;
        l=length(ans_zs(1,:));
        while i~=l-1
            if ans_zs(2,i)*ans_zs(2,i+1)<0
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
%                     fprintf('\n t = %.5f, f(t) = %.5f \n', c, f(c));
                    f_and_c=[f_and_c,[f1;c]];
                end
            end
            i=i+1;
        end
       %% 第二段
        ans_zs2=[];
        a=c1+1;
        b=cs0;
        %% 方程定义，模态方程按硬质海底条件展开(第二段)
        f = @(c) tan(2.*pi.*f1/c.*sqrt(c.^2/c1.^2-1).*h).*rou1/rou2.*c.^4/cs0.^4.*sqrt(1-c.^2/cp0.^2)/sqrt(c.^2/c1.^2-1)+(2-c.^2/cs0.^2).^2-4.*sqrt(1-c.^2/cs0.^2).*sqrt(1-c.^2/cp0.^2);
        for c=a:0.1:b
            ans_zs2=[ans_zs2,[c;f(c)]];
        end
        i=1;
        l=length(ans_zs2(1,:));
        while i~=l-1
            if ans_zs2(2,i).*ans_zs2(2,i+1)<0
                a=ans_zs2(1,i);
                b=ans_zs2(1,i+1);
                t=(a+b)./2;
                while abs(b-a)>1e-10
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
                end
            end
            i=i+1;
        end
    end
    % 这里调整共画几阶模态，需要调整数量时对应添加代码行
    f_and_c1 = [];
    f_and_c2 = [];
    
    for f2 = f_range
        [~,num] = find(f_and_c(1,:)==f2);
        f_and_c1 = [f_and_c1,[f_and_c(1,min(num));f_and_c(2,min(num))]];
        f_and_c2 = [f_and_c2,[f_and_c(1,min(num(2:end)));f_and_c(2,min(num(2:end)))]];
    end
    w_and_k1 = [2*pi*f_and_c1(1,:); (2*pi*f_and_c1(1,:))./f_and_c1(2,:)];
    w_and_k2 = [2*pi*f_and_c2(1,:); (2*pi*f_and_c2(1,:))./f_and_c2(2,:)];
    
    figure;
    scatter(f_and_c1(1,:),f_and_c1(2,:),5,'filled','b');
    hold on;
    scatter(f_and_c2(1,:),f_and_c2(2,:),5,'filled','b');
    scatter(f_and_c1(1,1:end-1),diff(w_and_k1(1,:))./diff(w_and_k1(2,:)),5,'filled','r');
    scatter(f_and_c2(1,1:end-1),diff(w_and_k2(1,:))./diff(w_and_k2(2,:)),5,'filled','r');
    hold off;
    
    title(['h=',num2str(h),' c_w=',num2str(c1),' \rho_w=',num2str(rou1),' c_p=',num2str(cp0),' c_s=',num2str(cs0),' \rho_e=',num2str(rou2),]);
    xlabel('f(Hz)');ylabel('c(m/s)')
end
