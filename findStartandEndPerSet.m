function [arclength,br,p_proj_sel,long_ip] = findStartandEndPerSet(points,radii,center,angle_eps)  
    br=0;p_proj_sel=[];long_ip=[];
    %点到圆的投影，用极坐标解决 
    [thetas,~]=cart2pol(points(:,1)-center(1),points(:,2)-center(2));
    %把极坐标转为直角坐标，因为基于密度的聚类用的直角坐标
    [x_on_pre,y_on_pre]=pol2cart(thetas,repmat(radii,1,length(thetas))');
    x_on_circle=x_on_pre+center(1);
    y_on_circle=y_on_pre+center(2);
    p_proj=[x_on_circle,y_on_circle];

    %epsilon;%可以和半径成正比，本质上就是固定度了。
    epsilon=2*radii*sin(angle_eps*pi/360); %15度
    minPts = 5;
    [idx, ~] = dbscan(p_proj, epsilon, minPts);%如果idx都是负数呢。
    if(max(idx)==-1)
        arclength=-1;
        return
    end
    %选出非噪声的类
    p_proj(:,3)=idx;
    p_proj_sel=p_proj(idx>0,:);
    %计算每个类的弧度范围，在极坐标中解决
    [thetas_sets,~]=cart2pol(p_proj_sel(:,1)-center(1),p_proj_sel(:,2)-center(2));
    thetas_sets(thetas_sets < 0) = thetas_sets(thetas_sets < 0) + 2*pi;
    thetas_sets(:,2)=p_proj_sel(:,3);
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
        %在直角坐标系中解决吧    %虽然有误差，但很小，把点按顺序链接起来的线段近似弧长
        order_thetas_sets=sortrows(thetas_sets,1);
        order_thetas_sets(:,2)=radii;
        [xx,yy]=pol2cart(order_thetas_sets(:,1),order_thetas_sets(:,2));
        xxyy=[xx,yy];
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
        range_a_set=arclength/radii;
        if abs(range_a_set-range_per_set(1,2))>0.02 %经过了0度 临界值是0.01722以上，但不能太大
            rows=length(order_thetas_sets);
            % 找到thetas变化最大的两个点，两个点的平均作为分割点
            diff=order_thetas_sets(:,1)-order_thetas_sets([2:rows,1],1);
            [~, minIndex] = min(diff);
            fenge=(order_thetas_sets(minIndex,1)+order_thetas_sets(minIndex+1,1))/2;
            largeset=thetas_sets(thetas_sets(:,1)>fenge,:);
            smallset=thetas_sets(thetas_sets(:,1)<fenge,:);
            range_breakset=2*pi-min(largeset(:,1))+max(smallset(:,1));
            arclength=range_breakset*radii;
            range_per_set(1,2)=range_breakset; %修正BR
        end
        long_ip=p_proj_sel;
    else %当存在一个以上的类
        if (sum(range_per_set(:,2))>2*pi-angle_eps*pi/180) %有一个类经过0度，那么范围和一定大于2pi-epsilon对应的弧度:15*pi/180=0.2618
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
            long_ip=p_proj_sel(p_proj_sel(:,3)==max_id,:);
            %选择弧长范围最大的类
            arclength=max([range_per_remainset(:,2)',range_breakset])*radii;
        else
            arclength=max(range_per_set(:,2))*radii; %没一个类经过0度，那么范围和一定小于2pi,直接找最大值即可
            [~,max_id]=max(range_per_set(:,2));
            long_ip=p_proj_sel(p_proj_sel(:,3)==max_id,:);
        end
    end
    br=max(range_per_set(:,2));
end
