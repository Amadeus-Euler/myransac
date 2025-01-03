function [r,p,ip_on_circle_sel,longest_ip_on_circle_sel]=myransac_plus2(data,sigma,itermax,max_r)
% ip_on_circle是内点在拟合圆上的投影并去除了外点，epsilon是20度对应的弦长，max_r是半径上限
r=0;p=0;%初始赋值
ip_on_circle_sel=[];
longest_ip_on_circle_sel=[];
a = data;
% RANSCA参数
% 迭代次数
iter = 0;
% 查看圆数据的大小
m=length(a);
if m<6
    fprintf('点数太少')
    r=0;p=0;
    return
end
% 误差参数

largestarclength = 0;

% 开始循环迭代
while iter<itermax

    % 随机挑选三个点，三个点不重复
    % 拟合圆最少需要三个点，拟合直线最少需要两个
    % ran为索引编号
    ran = randperm(m,3)';
    % b为索引得到的点
    b = a(ran,:);
    % 根据随机得到的三个点，计算圆的半径和圆心
    [r1,p1] = ThreePoint2Circle(b(1,1:2), b(2,1:2), b(3,1:2));
    if r1>max_r*1.1|| r1<0.01%设置上下限
        iter=iter+1;
        continue
    end
%     % 选择除了随机得到的三个点外的其他点
%     c = setdiff(a,b,"rows");
%   构成圆的点也算进内点
    c=a;
    % 计算每个点到圆心的距离dis
    dis = sqrt(sum((c(:,1:2)-p1).^2,2));
    % 计算 dis和拟合圆的误差
    res = abs(dis - r1);
    % 选择小于误差的点，进入到内点中
    d = c(res<sigma,:);
%         scatter(a(:,1),a(:,2),'MarkerEdgeColor','blue');
%         axis equal;
%         hold on
%         xfit = p1(1) + r1 * cos(theta);
%         yfit = p1(2) + r1 * sin(theta);
%         plot(xfit, yfit, '--', 'Color', 'red', 'LineWidth', 1.5);
%         scatter(d(:,1),d(:,2),'MarkerEdgeColor','red');
    if(isempty(d))
        iter=iter+1;
        continue
    end
    %点到圆的投影，用极坐标解决
    [thetas,~]=cart2pol(d(:,1)-p1(1),d(:,2)-p1(2));
    %把极坐标转为直角坐标，因为基于密度的聚类用的直角坐标
    [x_on_pre,y_on_pre]=pol2cart(thetas,repmat(r1,1,length(thetas))');
    x_on_circle=x_on_pre+p1(1);
    y_on_circle=y_on_pre+p1(2);
    d_proj=[x_on_circle,y_on_circle];

    %epsilon;%可以和半径成正比，本质上就是固定度了。但这个是弧长,后面再改改,r1*2*pi/18
    epsilon=2*r1*sin(15*pi/360); %15度
    minPts = 5;
    [idx, ~] = dbscan(d_proj, epsilon, minPts);%如果idx都是负数呢。
    if(max(idx)==-1)
        iter=iter+1;
        continue
    end
    %选出非噪声的类
%     d(:,3)=idx;
%     d_sel=d(idx>0,:);
    d_proj(:,3)=idx;
    d_proj_sel=d_proj(idx>0,:);

%         figure
%         scatter(x_on_circle,y_on_circle,'MarkerEdgeColor','red');
%         figure
%         scatter(d_proj(:,1), d_proj(:,2), 20, idx, 'filled');
    %计算每个类的弧度范围，在极坐标中解决
    [thetas_sets,~]=cart2pol(d_proj_sel(:,1)-p1(1),d_proj_sel(:,2)-p1(2));
    thetas_sets(thetas_sets < 0) = thetas_sets(thetas_sets < 0) + 2*pi;
    thetas_sets(:,2)=d_proj_sel(:,3);
    lengths=unique(thetas_sets(:, 2));
    range_per_set = zeros(length(lengths), 3);
    for j = 1:length(lengths)
        current_class = lengths(j);
        set_scores = thetas_sets(thetas_sets(:, 2) == current_class, 1);
        set_range = max(set_scores) - min(set_scores);
        set_sd=std(set_scores);
        range_per_set(j, :) = [current_class, set_range, set_sd];
    end
    %当只有一个类
    if (max(idx)==1)
        %仅一类且经过0度怎么判断？
        %arclength=range_per_set(:,2)*r1;
        %在直角坐标系中解决吧
        %虽然有误差，但很小，把点按顺序链接起来的线段近似弧长
        order_thetas_sets=sortrows(thetas_sets,1);
        order_thetas_sets(:,2)=r1;
        [xx,yy]=pol2cart(order_thetas_sets(:,1),order_thetas_sets(:,2));
        xxyy=[xx,yy];
        % scatter(xx, yy, 20, 'filled');
        %计算每个线段长度并求和
        xxyy_cuo=xxyy([length(xxyy),1:(length(xxyy)-1)],:);
        diff_mat=xxyy-xxyy_cuo;
        arclengths=zeros(length(xxyy),1);
        for k=1:length(xxyy)
            arclengths(k)=norm(diff_mat(k,:));
        end
        
        %有时候，整个类的弧长小于eplison,导致没有发现跳跃点，导致弧长大了两倍，删除最大的连线即可
        if sum(arclengths>epsilon)==0
            arclength=sum(arclengths)-max(arclengths);
        else
            arclengths(arclengths>epsilon)=0;  %大于epsilon的线段是跳跃点的连线，应去除。
            arclength=sum(arclengths);
        end
       
        range_a_set=arclength/r1;
%         if abs(range_a_set-range_per_set(1,2))>0.1 %经过了0度
%             break
%         end
        if abs(range_a_set-range_per_set(1,2))>0.02 %经过了0度 临界值是0.01722以上，但不能太大
            rows=length(order_thetas_sets);
            % 找到thetas变化最大的两个点，两个点的平均作为分割点
            diff=order_thetas_sets(:,1)-order_thetas_sets([2:rows,1],1);
            [~, minIndex] = min(diff);
            fenge=(order_thetas_sets(minIndex,1)+order_thetas_sets(minIndex+1,1))/2;
            largeset=thetas_sets(thetas_sets(:,1)>fenge,:);
            smallset=thetas_sets(thetas_sets(:,1)<fenge,:);
            range_breakset=2*pi-min(largeset(:,1))+max(smallset(:,1));
            arclength=range_breakset*r1;
        end
        long_ip=d_proj_sel;
        %当存在一个以上的类
    else
        if (sum(range_per_set(:,2))>2*pi-2*asin(epsilon / (2 * r1))) %有一个类经过0度，那么范围和一定大于2pi-epsilon对应的弧度:15*pi/180=0.2618
            maxrange=max(range_per_set(:,2));%范围最大的set经过了0度，找到范围最大的set
            maxrange_idx=range_per_set(range_per_set(:,2)==maxrange,1);%范围最大的set的idx为maxrange_idx
            remain_idx=range_per_set(range_per_set(:,2)<maxrange,1);%从剩余set中选择一个点作为分界点
            sets2selpoint=thetas_sets(ismember(thetas_sets(:,2),remain_idx));
            fenge=sets2selpoint(1);%分界点找到了
            sets2break=thetas_sets(thetas_sets(:,2)==maxrange_idx,:);
            %把经过了0度的类分为两类
            largeset=sets2break(sets2break(:,1)>fenge,:);
            smallset=sets2break(sets2break(:,1)<fenge,:);
            range_breakset=2*pi-min(largeset(:,1))+max(smallset(:,1));
            %剩余的类的类范围用原来的表去除经过了0度的类，找到
            range_per_remainset=range_per_set(range_per_set(:,1)~=maxrange_idx,:);
            range_per_set=[range_per_remainset; %这时第一列可能不是递增的
            [maxrange_idx,range_breakset,2]
            ];
            max_id=range_per_set(range_per_set(:,2)==max(range_per_set(:,2)),1);%选择最大范围的id
            long_ip=d_proj_sel(d_proj_sel(:,3)==max_id,:);
            %选择弧长范围最大的类
            arclength=max([range_per_remainset(:,2)',range_breakset])*r1;

        else
            arclength=max(range_per_set(:,2))*r1; %没一个类经过0度，那么范围和一定小于2pi,直接找最大值即可
            [~,max_id]=max(range_per_set(:,2));
            long_ip=d_proj_sel(d_proj_sel(:,3)==max_id,:);
        end
    end
%     scatter(d(:,1), d(:,2),10, 'MarkerEdgeColor','red');
%     hold on
%     scatter(long_ip(:,1), long_ip(:,2),30, 'MarkerEdgeColor','blue');


    % 存最优的圆
    if (arclength>=largestarclength)
        largestarclength = arclength;
        r = r1; p = p1; ip_on_circle_sel = d_proj_sel; longest_ip_on_circle_sel = long_ip;
    end
    iter = iter + 1;
end
% % 
%     scatter(ip_on_circle_sel(:,1), ip_on_circle_sel(:,2),10, 'MarkerEdgeColor','red');
%     hold on
%     scatter(longest_ip_on_circle_sel(:,1), longest_ip_on_circle_sel(:,2),30, 'MarkerEdgeColor','blue');





