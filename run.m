%% 10cm提取一次
% H=18.6;
alg='ols';
fenli=mat_sel;
theta = linspace(0, 2*pi, 100);
max_r=0.20;%25cm    

hs = fagen:0.1:H;
sigma=0.004;
ratios=zeros(length(hs),1);
ranges=zeros(length(hs),1);
ds = zeros(length(hs),1);
ps = zeros(length(hs),2);
% BRs = zeros(length(hs),1);
thick=0.02; %切片厚度
for i=1:length(hs)
    if mod(i-1,20) == 0
        figure
    end
    span = [hs(i)-thick/2,hs(i)+thick/2];
    sub_fenli = fenli(span(1)<fenli(:,3)&span(2)>fenli(:,3),:);
    crosss = sub_fenli(:,1:2);
    crosss = unique(crosss,'rows');
    if (isempty(crosss))
        continue
    end

    
    if strcmp(alg,'plus')
        [r,p,BR,ip_on_circle_sel,longest_ip_on_circle_sel] = myRanSac(crosss,sigma,300,max_r);
    end
    if strcmp(alg,'ran')
        [r,p,BR,ip_on_circle_sel,longest_ip_on_circle_sel] = RanSac(crosss,sigma,300,max_r);
    end
    if strcmp(alg,'ols')
        [r,p,BR,ip_on_circle_sel,longest_ip_on_circle_sel] = CircleFitByPratt3(crosss,sigma);
    end

    if(isempty(ip_on_circle_sel))
        continue
    end
    ip_on_circle_sel_t=array2table(ip_on_circle_sel,'VariableNames',{'x','y','set'});
    k=mod(i,20);
    if k==0
        k=20;
    end
    subplot(4,5,k);
    scatter(crosss(:,1),crosss(:,2),20,'filled');
    axis equal;
    hold on
    scatter(ip_on_circle_sel(:,1), ip_on_circle_sel(:,2),10, 'MarkerEdgeColor','red');
    scatter(longest_ip_on_circle_sel(:,1), longest_ip_on_circle_sel(:,2),30, 'MarkerEdgeColor','blue');
    if r==0
        fprintf([',序号',num2str(i),', 绝对高',num2str(hs(i)),'拟合失败\n']);
    end
    if r~=0
        xfit = p(1) + r * cos(theta);
        yfit = p(2) + r * sin(theta);
        plot(xfit, yfit, '--', 'Color', 'red', 'LineWidth', 1.5);
        title(['no:',num2str(i),', h=',num2str(round(hs(i),3)),', d=',num2str(round(r*200,3))]);
        hold off
        ds(i)=r*200;
        ps(i,:)=p;
    else
        title(['no:',num2str(i),', h=',num2str(hs(i)),', d=',num2str(r*200)]);
        hold off
    end
       %限制下一次直径的最大值
    max_r=median([ds(ds~=0);50])/200;
    ranges(i)=BR;
end
    hds=[hs',ds,ranges];
    hds_t=array2table(hds,"VariableNames",{'h','d','range'});



%% 最大弧度
hds_t_ex0=hds_t(hds_t.d>0,:);
hds_t_ex0_higher5=hds_t_ex0(hds_t_ex0.h>0.6,:);
heights=height(hds_t_ex0_higher5);
hds_t_ex0_higher5.save=zeros(height(hds_t_ex0_higher5),1);
interval=10;


for i=1:ceil(heights/interval)
    hds_sub=hds_t_ex0_higher5(i*interval-interval+1:min(heights,i*interval),:);  %最后一个可能没有5行
    sortedDf = sortrows(hds_sub,"range",'descend');
    choosed_h=sortedDf.h(1); %第一大范围的高度
    rowIndices = find(hds_t_ex0_higher5.h == choosed_h);
    hds_t_ex0_higher5.save(rowIndices)=1;
end
hds_t_ex0_higher5_sel=hds_t_ex0_higher5(hds_t_ex0_higher5.save==1,:);
hds_t_ex0_higher5_sel.save=[];
maxrange=[hds_t_ex0(hds_t_ex0.h<=0.6,:);hds_t_ex0_higher5_sel;{H,0,0}];%添加最高点
maxrange.save=ones(height(maxrange),1);

fitdata=maxrange(maxrange.save==1,:);
figure
scatter(fitdata.h,fitdata.d, [],fitdata.range);
lightBlue = [0.7, 0.85, 1]; % 浅蓝色 RGB
blue = [0, 0, 1];           % 蓝色 RGB
N = 256; % 色阶数
customColormap = [linspace(lightBlue(1), blue(1), N)', ...
                  linspace(lightBlue(2), blue(2), N)', ...
                  linspace(lightBlue(3), blue(3), N)'];
colormap(customColormap); % 设置自定义颜色映射

% 调整色阶范围
hold on
ft = fittype('smoothingspline');

opts = fitoptions(ft);
% opts.Weights = fitdata.ratio;
opts.SmoothingParam = 0.6; % 平滑参数 (惩罚系数) 范围为 [0(最为光滑), 1]
% 拟合数据
[fitresult, gof] = fit(fitdata.h, fitdata.d, ft, opts);
% 绘制结果
plot(fitresult);
ylim([0, 40]); 
title('最大弧度筛选三次样条回归');xlabel('X');ylabel('Y');
figure
scatter(hds_t_ex0.h,hds_t_ex0.d,[],hds_t_ex0.range);