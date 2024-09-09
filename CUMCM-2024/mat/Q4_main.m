clc, clear,close all
%% 程序主代码区 %%
%%%%%%%%%%%%%%%%%

% 获取最小调头半径, 考虑碰撞和速度限制 (实际是否能通过)
    t_step = 0.05;
    R_range_init = [0, 4.5];
    R_step1 = 0.5;
    R_step2 = 0.1;
    error_max_abs = 10^(-12);
    %R_t_min = 4.254673736180;
     R_t_min = 4.25467374;
    [R_t_min, ~] = GetMinRadius(t_step, R_range_init, R_step1, R_step2, error_max_abs);


% 根据此调头半径计算所需参数
    % 给定 R_turn_min
    global R_t; R_t = R_t_min;

    ChangeGlobalVaribles_Q4

% 由题目直接确定的参数
global d d_0 luoju b v_0 theta_max
    d = 1.65;    % 普通节点间距
    d_0 = 2.86;    % 龙头节点间距
    luoju = 1.7; % 螺距
    a = 0; % 阿基米德螺线的起始半径
    b = luoju / (2*pi); % 每圈的增长量 (与螺距有关)
    theta_max = 16 * 2 * pi;      % 初始点极坐标角度
    v_0 = 1;

if 0
% 依次求解全部所需全局变量
    global theta_t coor_M coor_H coor_E dist_HE l_M l_H t_turn ...
           l_HA coor_A coor_C R_small R_big lambda coor_HA kappa
        theta_t = R_t/b;
        func = @(x)b*sqrt(x.^2+1);
        t_turn = 1/v_0*integral(func, theta_t , 32*pi); % 从进入开始计时, 到恰好开始转弯需要的时间
        coor_M = b*theta_t*[cos(theta_t), sin(theta_t)];
        coor_H = -coor_M;
        coor_E = coor_H/3;
        dist_HE = R_t/3;
        l_M = - [cos(theta_t) - theta_t.*sin(theta_t),  sin(theta_t) + theta_t.*cos(theta_t)];
        l_H = - l_M;
        l_HA = - cross([0,0,1], [l_H, 0]); l_HA = l_HA(1:2); l_HA = l_HA/norm(l_HA);    % cross 前补了个负号
        func = @(t) norm(coor_E - (coor_H + t*l_HA)) - norm(coor_H - (coor_H + t*l_HA));
        t = fzero(func, 0);
        coor_A = coor_H + t*l_HA;
        coor_C = 3*coor_E - 2*coor_A;
        R_small = norm(coor_A - coor_H);
        R_big = 2*R_small;
        lambda = acos(1 - norm(coor_H-coor_E)^2/(2*R_small^2));
        coor_HA = coor_A - coor_H;
        kappa = pi - atan(coor_HA(2)/coor_HA(1)); 
        % 注意 kappa 的范围, 需要作判断 if(-coor_HA(1) < 0){ kappa --> kappa - pi}以纠正
         if ( - coor_HA(1) < 0)
             kappa = kappa - pi;
         end
    
    % 构建全局函数
    global P_C P_A
        % 圆 C 的参数方程 (返回直角坐标)
        P_C = @(beta) ( coor_C + R_big *cos(beta + kappa).*[1 0] - R_big *sin(beta + kappa).*[0 1] );
        % 圆 A 的参数方程 (返回直角坐标)
        P_A = @(beta) ( coor_A - R_small *cos((-beta) + kappa + lambda).*[1 0] + R_small *sin((-beta)  + kappa + lambda ).*[0 1]);
    
    % 用于计算距离的全局函数族
    global func_dist_21 func_dist_22 func_dist_32 func_dist_33 func_dist_42 func_dist_43 func_dist_44
        func_dist_21 = @(beta, theta, d) ... % 转直角后用距离公式, 并作映射 dist --> 2*d - dist 以满足 min 时 > d
            2*d - norm(P_C(beta) - b*theta.*[cos(theta), sin(theta)]);
        func_dist_22 = @(beta_1, beta_2) R_big * sqrt( 2*(1-cos(beta_1 - beta_2)) ); % 半径是 R_big
        func_dist_32 = @(beta_1, beta_2) norm(P_A(beta_1) - P_C(beta_2)); % 注意 beta_1 是 P_A
        func_dist_33 = @(beta_1, beta_2) R_small * sqrt( 2*(1-cos(beta_1 - beta_2)) ); % 半径是 R_small
        func_dist_42 = @(theta, beta) ...   % 这里给的 theta 转直角后再对称 (取负) 才是真实坐标
            norm( - b*theta.*[cos(theta), sin(theta)] - P_C(beta));
        func_dist_43 = @(theta, beta) ...   % 这里给的 theta 转直角后再对称才是真实坐标
            norm( - b*theta.*[cos(theta), sin(theta)] - P_A(beta));
        func_dist_44 = @(theta_1, theta_2) ...
            norm( b*theta_1.*[cos(theta_1), sin(theta_1)] - b*theta_2.*[cos(theta_2), sin(theta_2)]);
end
    
% 检查 R_turn_min 是否可行
    
    [logic_crush, logic_speederror, t_crushing] = CheckCrushAndSpeed_Q4([0, 40], 0.1);
    if (logic_crush == true) || (logic_speederror == true)
        disp(['R_t = ', num2str(R_t, '%.12f')]);
        error(">> error:  当前 R_t 发生碰撞或无法通过, 须重新设置 R_t 值")
    end

% 输出 Q4 最终结果
    t_array = -100:1:100;
    coordinates_matrix = zeros(224*2, length(t_array));
    speed_matrix = zeros(224, length(t_array));
    for j = 1:length(t_array)
        t = t_array(j);
        % 求解此 t 下的直角坐标
        angel_withflag_0 = GetHeadLocation(v_0, t);
        angle_withflag_all = GetAllPoints_Q4(angel_withflag_0, 10^(-16));
        coordinates_array = Angle2Coor(angle_withflag_all);
        speed_matrix(:, j) = GetAllSpeed_Q4(angle_withflag_all, v_0);
        % 将坐标格式转为 xlsx 中的模版格式
        for i = 1:2*224
            if mod(i, 2) == 1   % i 奇数, 对应 x 
                coordinates_matrix(i, j) = coordinates_array(ceil(i/2), 1);
            else                % i 偶数, 对应 y 
                coordinates_matrix(i, j) = coordinates_array(i/2, 2);
            end
        end
    end
    PrintResult_Q4(coordinates_matrix, speed_matrix)

% 可视化部分
    % 路径可视化
        DrawWholePath   
    % t = 100 结果可视化 
        t = 100;
        angle_withflag_0 = GetHeadLocation(v_0, t);
        angle_withflag_all = GetAllPoints_Q4(angle_withflag_0, 10^(-16));
        coor_all = Angle2Coor(angle_withflag_all);
        coor_Rectangles = GetRectangles_Q4(coor_all);
        DrawPoints_Q4(coor_all);
        DrawPointsAndRectangles_Q4(coor_all, coor_Rectangles);


%% 问题四函数区 %%
%%%%%%%%%%%%%%%%%


function ChangeGlobalVaribles_Q4
% 引入全局变量
global R_t d d_0 luoju b v_0 theta_max
% 依次求解全部所需全局变量
    global theta_t coor_M coor_H coor_E dist_HE l_M l_H t_turn ...
           l_HA coor_A coor_C R_small R_big lambda coor_HA kappa

        theta_t = R_t/b;
        func = @(x)b*sqrt(x.^2+1);
        t_turn = 1/v_0*integral(func, theta_t , 32*pi); % 从进入开始计时, 到恰好开始转弯需要的时间
        coor_M = b*theta_t*[cos(theta_t), sin(theta_t)];
        coor_H = -coor_M;
        coor_E = coor_H/3;
        dist_HE = R_t/3;
        l_M = - [cos(theta_t) - theta_t.*sin(theta_t),  sin(theta_t) + theta_t.*cos(theta_t)];
        l_H = - l_M;
        l_HA = - cross([0,0,1], [l_H, 0]); l_HA = l_HA(1:2); l_HA = l_HA/norm(l_HA);    % cross 前补了个负号
        func = @(t) norm(coor_E - (coor_H + t*l_HA)) - norm(coor_H - (coor_H + t*l_HA));
        t = fzero(func, 0);
        coor_A = coor_H + t*l_HA;
        coor_C = 3*coor_E - 2*coor_A;
        R_small = norm(coor_A - coor_H);
        R_big = 2*R_small;
        lambda = acos(1 - norm(coor_H-coor_E)^2/(2*R_small^2));
        coor_HA = coor_A - coor_H;
        kappa = pi - atan(coor_HA(2)/coor_HA(1)); 
        % 注意 kappa 的范围, 需要作判断 if(-coor_HA(1) < 0){ kappa --> kappa - pi}以纠正
         if ( - coor_HA(1) < 0)
             kappa = kappa - pi;
         end
    
    % 构建全局函数
    global P_C P_A
        % 圆 C 的参数方程 (返回直角坐标)
        P_C = @(beta) ( coor_C + R_big *cos(beta + kappa).*[1 0] - R_big *sin(beta + kappa).*[0 1] );
        % 圆 A 的参数方程 (返回直角坐标)
        P_A = @(beta) ( coor_A - R_small *cos((-beta) + kappa + lambda).*[1 0] + R_small *sin((-beta)  + kappa + lambda ).*[0 1]);
    
    % 用于计算距离的全局函数族
    global func_dist_21 func_dist_22 func_dist_32 func_dist_33 func_dist_42 func_dist_43 func_dist_44
        func_dist_21 = @(beta, theta, d) ... % 转直角后用距离公式, 并作映射 dist --> 2*d - dist 以满足 min 时 > d
            2*d - norm(P_C(beta) - b*theta.*[cos(theta), sin(theta)]);
        func_dist_22 = @(beta_1, beta_2) R_big * sqrt( 2*(1-cos(beta_1 - beta_2)) ); % 半径是 R_big
        func_dist_32 = @(beta_1, beta_2) norm(P_A(beta_1) - P_C(beta_2)); % 注意 beta_1 是 P_A
        func_dist_33 = @(beta_1, beta_2) R_small * sqrt( 2*(1-cos(beta_1 - beta_2)) ); % 半径是 R_small
        func_dist_42 = @(theta, beta) ...   % 这里给的 theta 转直角后再对称 (取负) 才是真实坐标
            norm( - b*theta.*[cos(theta), sin(theta)] - P_C(beta));
        func_dist_43 = @(theta, beta) ...   % 这里给的 theta 转直角后再对称才是真实坐标
            norm( - b*theta.*[cos(theta), sin(theta)] - P_A(beta));
        func_dist_44 = @(theta_1, theta_2) ...
            norm( b*theta_1.*[cos(theta_1), sin(theta_1)] - b*theta_2.*[cos(theta_2), sin(theta_2)]);
end


function [R_t_min, t_crushing] = GetMinRadius(t_step, R_range_init, R_step1, R_step2, error_max_abs)
global R_t v_0 b
    disp('函数 GetMinRadius 运行时会修改几乎所有全局变量')
    disp('需在运行后对全局变量重新赋值!')
    R_range = R_range_init
    
    start = tic;

    % 第一层步长: R_step1
    R_array = R_range(1):R_step1:R_range(2);
        for i = length(R_array):(-1):1
            R_t = R_array(i);
            %func = @(x)b*sqrt(x.^2+1);
            %t_turn_now = 1/v_0*integral(func, R_t/b , 32*pi); % 从进入开始计时, 到恰好开始转弯需要的时间
            %[logic_crush, logic_speederror, t_crushing] = CheckCrushAndSpeed_Q4([0, t_turn_now/3], t_step);
            [logic_crush, logic_speederror, t_crushing] = CheckCrushAndSpeed_Q4([0, 40], t_step);
            disp(['crush = ', num2str(logic_crush), ', speederror = ', num2str(logic_speederror), ', R = ', num2str(R_t), ', t_crushing = ', num2str(t_crushing)]);
            if (logic_crush == true) || (logic_speederror == true)
                R_range = [R_array(i), R_array(i+1)]
                break
            end
        end

    % 第二层步长: R_step2
    R_array = R_range(1):R_step2:R_range(2);
        for i = length(R_array):(-1):1
            R_t = R_array(i);
            %func = @(x)b*sqrt(x.^2+1);
            %t_turn_now = 1/v_0*integral(func, R_t/b , 32*pi); % 从进入开始计时, 到恰好开始转弯需要的时间
            %[logic_crush, logic_speederror, t_crushing] = CheckCrushAndSpeed_Q4([0, t_turn_now/3], t_step);
            [logic_crush, logic_speederror, t_crushing] = CheckCrushAndSpeed_Q4([0, 40], t_step);
            disp(['crush = ', num2str(logic_crush), ', speederror = ', num2str(logic_speederror), ', R = ', num2str(R_t), ', t_crushing = ', num2str(t_crushing)]);
            if (logic_crush == true) || (logic_speederror == true)
                R_range = [R_array(i), R_array(i+1)]
                break
            end
        end

    % 牛顿迭代法
    min = R_range(1);
    max = R_range(2);
    for i = 1:71     % ceil( (log(4*pi)+20*log(10))/log(2) ) = 71, 这使得 4*pi/2^n < 10^(-20)
        error_max = max - min;
        temp = 0.5 * (max + min);
        R_t = temp;
        %func = @(x)b*sqrt(x.^2+1);
        %t_turn_now = 1/v_0*integral(func, R_t/b , 32*pi); % 从进入开始计时, 到恰好开始转弯需要的时间
        %[logic_crush, logic_speederror, t_crushing] = CheckCrushAndSpeed_Q4([0, t_turn_now/3], t_step);
        [temp_logic_crush, temp_logic_speederror, temp_t] = CheckCrushAndSpeed_Q4([0, 40], t_step);
        if (temp_logic_crush == true) || (temp_logic_speederror == true)
            min = temp
            t_crushing = temp_t;
        else
            max = temp
        end
        disp(['crush = ', num2str(temp_logic_crush), ', speederror = ', num2str(temp_logic_speederror), ', R = ', num2str(R_t), ', t_crushing = ', num2str(t_crushing)]);
        if error_max <= error_max_abs
            break
        end
    end

    % 记录结果
    R_t_min = 0.5*(min + max);
    time = toc(start);
    output_cell{1, 1} = num2str(time, '%.6f');  
    output_cell{1, 2} = num2str(R_t_min, '%.12f'); 
    output_cell{1, 3} = num2str(t_crushing, '%.12f');
    output_cell{1, 4} = num2str(t_step, '%.6f');
    output_cell{1, 4} = num2str(R_range_init, '%.6f');
    output_cell{1, 5} = num2str(R_step1, '%.6f');
    output_cell{1, 6} = num2str(R_step2, '%.6f');
    output_cell{1, 7} = num2str(error_max_abs, '%.20f');
    disp(['time = ', num2str(time, '%.6f')])
    disp(['R_t_min = ', num2str(R_t_min, '%.12f')])
    disp(['t_crushing = ', num2str(t_crushing, '%.6f')])
    disp(['t_step = ', num2str(t_step, '%.4f')])
    disp(['R_range_init = [', num2str(R_range_init(1), '%.3f'), ' ,', num2str(R_range_init(2), '%.3f'), ']'])
    disp(['R_step1 = ', num2str(R_step1, '%.4f')])
    disp(['R_step2 = ', num2str(R_step2, '%.4f')])
    disp(['R_error_max_abs = ', num2str(error_max_abs, '%.16f')])
    writecell(output_cell, 'MinRadius_Result.xlsx', 'WriteMode', 'append');
end    

function [logic_crush, logic_speederror, t_crushing] = CheckCrushAndSpeed_Q4(t_range, t_step)

    % 外层函数修改了 R_turn
    global R_t;
    
    % 题目确定的参数
    % 参数设置
    global d d_0 luoju b v_0 theta_max
        d = 1.65;    % 普通节点间距
        d_0 = 2.86;    % 龙头节点间距
        luoju = 1.7; % 螺距
        a = 0; % 阿基米德螺线的起始半径
        b = luoju / (2*pi); % 每圈的增长量 (与螺距有关)
        theta_max = 16 * 2 * pi;      % 初始点极坐标角度
        v_0 = 1;
    
    % 依次求解全部所需全局变量
    global theta_t coor_M coor_H coor_E dist_HE l_M l_H t_turn ...
           l_HA coor_A coor_C R_small R_big lambda coor_HA kappa
        theta_t = R_t/b;
        func = @(x)b*sqrt(x.^2+1);
        t_turn = 1/v_0*integral(func, theta_t , 32*pi); % 从进入开始计时, 到恰好开始转弯需要的时间
        coor_M = b*theta_t*[cos(theta_t), sin(theta_t)];
        coor_H = -coor_M;
        coor_E = coor_H/3;
        dist_HE = R_t/3;
        l_M = - [cos(theta_t) - theta_t.*sin(theta_t),  sin(theta_t) + theta_t.*cos(theta_t)];
        l_H = - l_M;
        l_HA = - cross([0,0,1], [l_H, 0]); l_HA = l_HA(1:2); l_HA = l_HA/norm(l_HA);    % cross 前补了个负号
        func = @(t) norm(coor_E - (coor_H + t*l_HA)) - norm(coor_H - (coor_H + t*l_HA));
        h = fzero(func, 0);
        coor_A = coor_H + h*l_HA;
        coor_C = 3*coor_E - 2*coor_A;
        R_small = norm(coor_A - coor_H);
        R_big = 2*R_small;
        lambda = acos(1 - norm(coor_H-coor_E)^2/(2*R_small^2));
        coor_HA = coor_A - coor_H;
        kappa = pi - atan(coor_HA(2)/coor_HA(1)); 
        % 注意 kappa 的范围, 需要作判断 if(-coor_HA(1) < 0){ kappa --> kappa - pi}以纠正
         if ( - coor_HA(1) < 0)
             kappa = kappa - pi;
         end
    
    % 构建全局函数
    global P_C P_A
        % 圆 C 的参数方程 (返回直角坐标)
        P_C = @(beta) ( coor_C + R_big *cos(beta + kappa).*[1 0] - R_big *sin(beta + kappa).*[0 1] );
        % 圆 A 的参数方程 (返回直角坐标)
        P_A = @(beta) ( coor_A - R_small *cos((-beta) + kappa + lambda).*[1 0] + R_small *sin((-beta)  + kappa + lambda ).*[0 1]);
    
    % 用于计算距离的全局函数族
    global func_dist_21 func_dist_22 func_dist_32 func_dist_33 func_dist_42 func_dist_43 func_dist_44
        func_dist_21 = @(beta, theta, d) ... % 转直角后用距离公式, 并作映射 dist --> 2*d - dist 以满足 min 时 > d
            2*d - norm(P_C(beta) - b*theta.*[cos(theta), sin(theta)]);
        func_dist_22 = @(beta_1, beta_2) R_big * sqrt( 2*(1-cos(beta_1 - beta_2)) ); % 半径是 R_big
        func_dist_32 = @(beta_1, beta_2) norm(P_A(beta_1) - P_C(beta_2)); % 注意 beta_1 是 P_A
        func_dist_33 = @(beta_1, beta_2) R_small * sqrt( 2*(1-cos(beta_1 - beta_2)) ); % 半径是 R_small
        func_dist_42 = @(theta, beta) ...   % 这里给的 theta 转直角后再对称 (取负) 才是真实坐标
            norm( - b*theta.*[cos(theta), sin(theta)] - P_C(beta));
        func_dist_43 = @(theta, beta) ...   % 这里给的 theta 转直角后再对称才是真实坐标
            norm( - b*theta.*[cos(theta), sin(theta)] - P_A(beta));
        func_dist_44 = @(theta_1, theta_2) ...
            norm( b*theta_1.*[cos(theta_1), sin(theta_1)] - b*theta_2.*[cos(theta_2), sin(theta_2)]);

    
    for t = t_range(1):t_step:t_range(2)
        angle_withflag_0 = GetHeadLocation(v_0, t);
        angle_withflag_all = GetAllPoints_Q4(angle_withflag_0, 10^(-16));
        coor_all = Angle2Coor(angle_withflag_all);
        logic_crush = IsCrushed_Q4(coor_all);
        if logic_crush == true
            t_crushing = t;
            logic_speederror = ~isempty(GetAllSpeed_Q4(angle_withflag_all, v_0)<0);
            return
        end
        speed_all = GetAllSpeed_Q4(angle_withflag_all, v_0);
        if ~isempty(find(speed_all<0, 1))
            t_crushing = t;
            logic_speederror = true;
            return
        end
    end

    % no crush or speed error
    logic_speederror = false;
    t_crushing = -1;
    return
end

function logic = IsCrushed_Q4(coor_all)
    % 数据准备
    d_special = sqrt(3.135^2 + 0.15^2) + sqrt(0.275^2 + 0.15^2);  % 龙头 3 号点最大判断间距
    d_common = sqrt(1.925^2 + 0.15^2) + sqrt(0.275^2 + 0.15^2);  % 2, 3 号点最大判断间距
    coor_Rectangles = GetRectangles_Q4(coor_all);

    % 与问题二不同，这里需要判断全部方框的四个点是否 legal

    % 判断龙头 (第一个板凳)
        point_four(1:4, 1:2) = coor_Rectangles(1:4, 1, :);

        % 挑选可能的节点
        distance_array = sqrt( sum( ( coor_all - coor_all(1,:) ).^2 , 2) );
        Logic = (distance_array < d_special) ;
        Logic(1:2) = 0; % 忽略第 1(本身),2(之后)
        Logic(224) = 0; % 忽略最后一个节点 (无板凳)
        for i = find(Logic)'    % 这里必须转置使得右值为行向量
            % 新的坐标原点
            new_Origin(1, 1:2) = coor_Rectangles(1, i, :);
            % 依次检查 4 个点
            for k = 1:4
                new_cor = (point_four(k, :) - new_Origin ) / [
                    coor_Rectangles(2, i, 1) - coor_Rectangles(1, i, 1), coor_Rectangles(2, i, 2) - coor_Rectangles(1, i, 2);
                    coor_Rectangles(4, i, 1) - coor_Rectangles(1, i, 1), coor_Rectangles(4, i, 2) - coor_Rectangles(1, i, 2);
                    ];   % 坐标系转换
                %disp(['龙头2号点 new_cor:', num2str(new_cor)])
                if (0<new_cor(1)&&new_cor(1)<1) && (0<new_cor(2)&&new_cor(2)<1)
                    logic = true;
                    %disp(['龙头: i = ',num2str(i),', k = ',num2str(k)])
                    return
                end
            end
        end

    % 判断其它板凳
        for j = 2:223
            point_four(1:4, 1:2) = coor_Rectangles(1:4, j, :);
    
            % 挑选可能的节点
            distance_array = sqrt( sum( ( coor_all - coor_all(1,:) ).^2 , 2) );
            Logic = (distance_array < d_common) ;
            Logic([j-1, j, j+1]) = 0; % 忽略第 j-1 (之前), j (本身), j+1(之后) 号板凳
            Logic(224) = 0; % 忽略最后一个节点 (无板凳)
            for i = find(Logic)'    % 这里必须转置使得右值为行向量
                % 新的坐标原点
                new_Origin(1, 1:2) = coor_Rectangles(1, i, :);
                % 依次检查 4 个点
                for k = 1:4
                    new_cor = (point_four(k, :) - new_Origin ) / [
                        coor_Rectangles(2, i, 1) - coor_Rectangles(1, i, 1), coor_Rectangles(2, i, 2) - coor_Rectangles(1, i, 2);
                        coor_Rectangles(4, i, 1) - coor_Rectangles(1, i, 1), coor_Rectangles(4, i, 2) - coor_Rectangles(1, i, 2);
                        ];   % 坐标系转换
                    if (0<new_cor(1)&&new_cor(1)<1) && (0<new_cor(2)&&new_cor(2)<1)
                        logic = true;
                        %disp(['j = ',num2str(j), ': i = ',num2str(i),', k = ',num2str(k)])
                        return
                    end
                end
            end
        end
        
    % 检查通过
    logic = false;
    return
end



function DrawPointsAndRectangles_Q4(coor_all, coor_Rectangles)
% 作线和点
    DrawPoints_Q4(coor_all);
% 作方框 
    coor_Rectangles(5, :, :) = coor_Rectangles(1, :, :);    % plot 围成闭合曲线
    hold on
    for i = 1:223
        plot(coor_Rectangles(:, i, 1), coor_Rectangles(:, i, 2), 'LineWidth', 0.3, 'Color', [1 0 1]);
    end
% 收尾
    hold off
end

function DrawWholePath
global theta_t b lambda P_C P_A R_t
    % 生成全路径曲线
    theta_array = theta_t:0.1:32*pi;
    beta_array = 0:0.05:lambda;
    coor_1_array = b*theta_array'.*[cos(theta_array'), sin(theta_array')];
    coor_4_array = - coor_1_array;
    coor_2_array = zeros(length(beta_array), 2);
    coor_3_array = zeros(length(beta_array), 2);
    for i = 1:length(beta_array)
        coor_2_array(i, :) = P_C(beta_array(i));
        coor_3_array(i, :)  = P_A(beta_array(i));
    end
    % 生成调头区域
    coor_5 = 4.5*[cos(0:0.02:2*pi); sin(0:0.02:2*pi)]'; % 给定调头区域
    coor_6 = R_t*[cos(0:0.02:2*pi); sin(0:0.02:2*pi)]'; % 实际调头区域

    % 作图
        figure('Color', [1 1 1])
        % 作出全路径曲线
        stc.line1 = plot(coor_1_array(:, 1), coor_1_array(:, 2));
        hold on
        stc.line23 = plot([coor_2_array(:, 1); coor_3_array(:, 1)], [coor_2_array(:, 2); coor_3_array(:, 2)]);
        stc.line4 = plot(coor_4_array(:, 1), coor_4_array(:, 2));
        stc.line5 = plot(coor_5(:, 1), coor_5(:, 2), 'black--');
        stc.line6 = plot(coor_6(:, 1), coor_6(:, 2), 'r--');
        
    % 设置样式
        % 坐标轴
            stc.fig = gcf;
            axis equal  
            stc.axes = gca;
            stc.axes.FontName = "Times New Roman"; % 全局 FontName
            stc.axes.Box = 'on';  
            stc.axes.FontSize = 14;
            xline(0, 'LineWidth', 0.3, 'Color',  [0.7, 0.7, 0.7]);
            yline(0, 'LineWidth', 0.3, 'Color',  [0.7, 0.7, 0.7]);
            stc.label.x = xlabel(stc.axes, '$x$', 'Interpreter', 'latex', 'FontSize', 15);
            stc.label.y = ylabel(stc.axes, '$y$', 'Interpreter', 'latex', 'FontSize', 15);
            

        % 标题
            stc.axes.Title.String = '';
            stc.axes.Title.FontSize = 17;
            stc.axes.Title.FontWeight = 'bold';
        
        % 线的样式
            stc.line1.LineWidth = 0.8;
            stc.line23.LineWidth = 0.8;
            stc.line4.LineWidth = 0.8;
            stc.line5.LineWidth = 0.6;
            stc.line6.LineWidth = 0.6; 
            stc.line1.Color = [0 1 0];   % 绿色
            stc.line23.Color = [1 0 0];   % 红色
            stc.line4.Color = [0 0 1];   % 蓝色
    
    % 收尾
        hold off
end


function Rectangle_Points = GetRectangles_Q4(coor_all)
    Vec_X = diff(coor_all);
    Vec_X = Vec_X ./ sqrt(sum(Vec_X.^2, 2));     % 单位化
    Vec_N = [ -Vec_X(: ,2), Vec_X(: ,1) ];  % 法向量


    % 计算矩形坐标
    Rectangle_P1 = coor_all(1:223, :) - Vec_X*0.275 + Vec_N*0.15;    % 注意是 1:223
    Rectangle_P2 = coor_all(1:223, :) - Vec_X*0.275 - Vec_N*0.15;
    Rectangle_P3 = coor_all(2:224, :) + Vec_X*0.275 - Vec_N*0.15;    % 注意是 2:224
    Rectangle_P4 = coor_all(2:224, :) + Vec_X*0.275 + Vec_N*0.15;

    Rectangle_Points = zeros(4, 223, 2);
    Rectangle_Points(1, :, :) = Rectangle_P1;
    Rectangle_Points(2, :, :) = Rectangle_P2;
    Rectangle_Points(3, :, :) = Rectangle_P3;
    Rectangle_Points(4, :, :) = Rectangle_P4;
end


function PrintResult_Q4(coor_matrix, speed_matrix)   
    % coor_matrix: (224*2)*(201) 矩阵
    % speed_matrix: 224*(201) 矩阵

    % 遍历矩阵，将每个元素格式化为保留 6 位小数的字符串  
    coor_xlsx = cell(size(coor_matrix));
    speed_xlsx = cell(size(speed_matrix));
    for i = 1:size(coor_matrix, 1)  
        for j = 1:size(coor_matrix, 2)  
            coor_xlsx{i, j} = num2str(coor_matrix(i, j), '%.6f');  
        end  
    end  
    for i = 1:size(speed_matrix, 1)  
        for j = 1:size(speed_matrix, 2)  
            speed_xlsx{i, j} = num2str(speed_matrix(i, j), '%.6f');  
        end  
    end  
    writecell(coor_xlsx, 'Q4_Result.xlsx', 'Sheet', 'Sheet1'); % 输出位置
    writecell(speed_xlsx, 'Q4_Result.xlsx', 'Sheet', 'Sheet2'); % 输出速度
end



function Speed_all = GetAllSpeed_Q4(angle_withflag_all, v_0)
% 引入全局变量
global P_C P_A coor_C coor_A
    Speed_all = zeros(224 ,1);
    Speed_all(1) = v_0;

    % 计算 Vec_X
    Coordinates = Angle2Coor(angle_withflag_all);
    Vec_X = diff(Coordinates);

    % 计算 Vec_Tao
    Vec_Tao = zeros(224, 2);
        % flag 1 上的点
        for i = find(angle_withflag_all(:, 2) == 1)'
            theta = angle_withflag_all(i, 1);
             Vec_Tao(i, :) = - [cos(theta) - theta*sin(theta), sin(theta) + theta*cos(theta) ];
        end
        % flag 2 上的点
        for i = find(angle_withflag_all(:, 2) == 2)'
            vec_r = [P_C(angle_withflag_all(i, 1)) - coor_C, 0];
            vec_tao = cross(vec_r, [0 0 1]);
            Vec_Tao(i, :) = vec_tao(1:2);
        end
        % flag 3 上的点
        for i = find(angle_withflag_all(:, 2) == 3)'
            vec_r = [P_A(angle_withflag_all(i, 1)) - coor_A, 0];
            vec_tao = cross([0 0 1], vec_r);    % [0 0 1] 在前
            Vec_Tao(i, :) = vec_tao(1:2);
        end
        % flag 4 上的点
        for i = find(angle_withflag_all(:, 2) == 4)'
            theta = angle_withflag_all(i, 1);
            Vec_Tao(i, :) = - [cos(theta) - theta*sin(theta), sin(theta) + theta*cos(theta) ];
        end

        Vec_Tao = Vec_Tao ./sqrt(sum( (Vec_Tao).^2 ,2));
    % 计算速度
    for i = 1:223
        % 速度正负
        mp = (   (sum((Vec_Tao(i, :).*Vec_X(i, :))) * sum((Vec_Tao(i+1, :).*Vec_X(i, :))))   >= 0    )*2 + -1; 
        Speed_all(i+1) = mp * Speed_all(i) * ( norm(Vec_Tao(i+1, :))*abs(sum(Vec_Tao(i, :).*Vec_X(i, :)))  ) / ( norm(Vec_Tao(i, :))*abs(sum(Vec_Tao(i+1, :).*Vec_X(i, :))) );
    end
end


function stc = DrawPoints_Q4(coordinates_array)
% 引入全局变量
global theta_t b lambda P_C P_A R_t

    % 生成全路径曲线
    theta_array = theta_t:0.1:32*pi;
    beta_array = 0:0.005:lambda;
    coor_1_array = b*theta_array'.*[cos(theta_array'), sin(theta_array')];
    coor_4_array = - coor_1_array;
    coor_2_array = zeros(length(beta_array), 2);
    coor_3_array = zeros(length(beta_array), 2);
    for i = 1:length(beta_array)
        coor_2_array(i, :) = P_C(beta_array(i));
        coor_3_array(i, :)  = P_A(beta_array(i));
    end
    % 生成调头区域
    coor_5 = 4.5*[cos(0:0.02:2*pi); sin(0:0.02:2*pi)]';
    coor_6 = R_t*[cos(0:0.02:2*pi); sin(0:0.02:2*pi)]'; % 实际调头区域
        

    % 作图
        figure('Color', [1 1 1])
        % 作出全路径曲线
        stc.line1 = plot(coor_1_array(:, 1), coor_1_array(:, 2));
        hold on
        stc.line23 = plot([coor_2_array(:, 1); coor_3_array(:, 1)], [coor_2_array(:, 2); coor_3_array(:, 2)]);
        stc.line4 = plot(coor_4_array(:, 1), coor_4_array(:, 2));
        stc.line5 = plot(coor_5(:, 1), coor_5(:, 2), 'black--');
        stc.line6 = plot(coor_6(:, 1), coor_6(:, 2), 'r--');
        % 作出所有节点
        stc.point_0 = scatter(coordinates_array(1, 1), coordinates_array(1, 2), 150, '.');
        stc.point_0.CData = [1 0 1];   % 粉色
        stc.points_rest = scatter(coordinates_array(2:end, 1), coordinates_array(2:end, 2), 50, 'black.');

    % 设置样式
        % 坐标轴
            stc.fig = gcf;
            axis equal  
            stc.axes = gca;
            stc.axes.FontName = "Times New Roman"; % 全局 FontName
            stc.axes.Box = 'on';  
            stc.axes.FontSize = 14;
            xline(0, 'LineWidth', 0.3, 'Color',  [0.7, 0.7, 0.7]);
            yline(0, 'LineWidth', 0.3, 'Color',  [0.7, 0.7, 0.7]);
            stc.label.x = xlabel(stc.axes, '$x$', 'Interpreter', 'latex', 'FontSize', 15);
            stc.label.y = ylabel(stc.axes, '$y$', 'Interpreter', 'latex', 'FontSize', 15);
            

        % 标题
            stc.axes.Title.String = '';
            stc.axes.Title.FontSize = 17;
            stc.axes.Title.FontWeight = 'bold';
        
        % 线的样式
            stc.line1.LineWidth = 0.8;
            stc.line23.LineWidth = 0.8;
            stc.line4.LineWidth = 0.8;
            stc.line5.LineWidth = 0.6;
            stc.line1.Color = [0 1 0];   % 绿色
            stc.line23.Color = [1 0 0];   % 红色
            stc.line4.Color = [0 0 1];   % 蓝色
    
    % 收尾
        hold off
end


function coordinates_array = Angle2Coor(angle_withflag_array)
% 将 angle_withflag_array (n*2 矩阵, 也即一列向量) 转为 直角坐标 coordinates (n*2 矩阵, 也即一列坐标向量)
% 引入全局变量
global P_C P_A b
    coordinates_array = zeros(size(angle_withflag_array));
    for i = 1 : size(angle_withflag_array, 1)
        switch angle_withflag_array(i, 2)
            case 1
                coordinates_array(i, :) = b*angle_withflag_array(i, 1) ...
                                           *[cos(angle_withflag_array(i, 1)), sin(angle_withflag_array(i, 1))];
            case 2
                coordinates_array(i, :) = P_C(angle_withflag_array(i, 1));
            case 3
                coordinates_array(i, :) = P_A(angle_withflag_array(i, 1));
            case 4
                coordinates_array(i, :) = - b*angle_withflag_array(i, 1) ...
                                            *[cos(angle_withflag_array(i, 1)), sin(angle_withflag_array(i, 1))];
        end
    end
end


function angle_withflag_all = GetAllPoints_Q4(angle_withflag_0, error_max_abs)
% 引入全局变量
global d d_0
    angle_withflag_all = zeros(224, 2);
    angle_withflag_all(1, :) = angle_withflag_0;
    % 龙头板凳长不同, 单独计算 第二节点
    angle_withflag_all(2, :) = GetNextPoint_newton_Q4(angle_withflag_all(1, :), d_0, error_max_abs);
    % 其它节点
    for i = 3:224
        angle_withflag_all(i, :) = GetNextPoint_newton_Q4(angle_withflag_all(i-1, :), d, error_max_abs);
    end
end


function angel_withflag = GetHeadLocation(v_0, t)
    % 引入全局变量
    global lambda R_big R_small b t_turn
    % 判断龙头 flag 标志位, 并给出对应 angle
    if t < 0
        angel_withflag(2) = 1;
        angel_withflag(1) = 16*2*pi - func_find_start(b, 16*2*pi, v_0*(t + t_turn));
    elseif ( 0 <= t ) && ( t < lambda*R_big/v_0 )
        angel_withflag(2) = 2;
        angel_withflag(1) = v_0*t/R_big;
    elseif ( lambda*R_big/v_0 <= t ) && ( t < lambda*(R_big+R_small)/v_0 )
        angel_withflag(2) = 3;
        angel_withflag(1) = v_0*(t - lambda*R_big/v_0) / R_small;
    elseif (t >= lambda*(R_big+R_small)/v_0)
        angel_withflag(2) = 4;
        angel_withflag(1) = 16*2*pi - func_find_start(b, 16*2*pi, v_0*(t_turn + lambda*(R_big + R_small)/v_0 - t));
    end
    return
end


function angle_withflag_next = GetNextPoint_newton_Q4(angle_withflag, d, error_max_abs)
    % 引入全局变量
    global b func_dist_21 func_dist_22 func_dist_32 func_dist_33 func_dist_42 func_dist_43 func_dist_44 theta_t lambda
    % 逻辑判断并求解
    switch angle_withflag(2)
        case 1
            angle_withflag_next = [GetNextPoint_newton(angle_withflag(1), b, d, error_max_abs), 1];
        case 2
            dist = func_dist_22(angle_withflag(1), 0);  % 计算距离
            if dist > d     % NextPoint 仍在 flag 2, 使用函数 func_dist_22
                % 其中 angle_1 定死, angle_2_range 用于遍历
                angle = GetPriciseAngle(angle_withflag(1), [0, angle_withflag(1)], d, error_max_abs, func_dist_22);
                angle_withflag_next = [angle, 2];
            else    % NextPoint 在 flag 1, 使用函数 func_dist_21
                angle = GetPriciseAngle(angle_withflag(1), [theta_t, theta_t + pi], d, error_max_abs, @(x,y) func_dist_21(x,y,d));
                angle_withflag_next = [angle, 1];
            end
        case 3
            dist = func_dist_33(angle_withflag(1), 0);  % 计算距离
            if dist > d     % NextPoint 仍在 flag 3, 使用函数 func_dist_33
                % 其中 angle_1 定死, angle_2_range 用于遍历
                angle = GetPriciseAngle(angle_withflag(1), [0, angle_withflag(1)], d, error_max_abs, func_dist_33);
                angle_withflag_next = [angle, 3];
            else    % NextPoint 在 flag 2, 使用函数 func_dist_32
                angle = GetPriciseAngle(angle_withflag(1), [0, lambda], d, error_max_abs, func_dist_32);
                angle_withflag_next = [angle, 2];
            end
        case 4
            dist = func_dist_44(angle_withflag(1), theta_t);  % 计算距离
            if (dist > d) || (angle_withflag(1) > theta_t + pi*2/3)     % NextPoint 仍在 flag 4, 使用函数 func_dist_44
                % 其中 angle_1 定死, angle_2_range 用于遍历
                % 这里 angle_range 选取不当会导致结果出错
                range_min = angle_withflag(1) - pi/3;
                if func_dist_44(angle_withflag(1), range_min) < d     % range_min 太大, 还需缩小
                    range_min = range_min - pi/3;
                    if func_dist_44(angle_withflag(1), range_min) < d
                        error("func_dist_44(angle_withflag(1), angle_withflag(1) - 2*pi/3) < d, range 传入错误!");
                    end
                end
                angle = GetPriciseAngle(angle_withflag(1), [range_min, angle_withflag(1)], d, error_max_abs, func_dist_44);
                angle_withflag_next = [angle, 4];
            else    % NextPoint 在 flag 3 或 flag 2
                if func_dist_43(angle_withflag(1), 0) < d
                    % NextPoint 在 flag 2, 使用函数 func_dist_42
                    angle = GetPriciseAngle(angle_withflag(1), [0, lambda], d, error_max_abs, func_dist_42);
                    angle_withflag_next = [angle, 2];
                else
                    % NextPoint 在 flag 3, 使用函数 func_dist_43
                    angle = GetPriciseAngle(angle_withflag(1), [0, lambda], d, error_max_abs, func_dist_43);
                    angle_withflag_next = [angle, 3];
                end
            end
    end
    return
end

function angel = GetPriciseAngle(angel_1, angle_2_range, d, error_max_abs, func_dist)
% func_dist(angle_1, angle_2, b), 其中 angle_1 定死, angle_2 用于遍历
% 注：此函数已完成去耦分离
% 函数要求 range 的 min 角度对应距离 > d, max 角度对应距离 < d
    % 第 1 层
    angle_array = linspace(angle_2_range(1), angle_2_range(2), 100);
    for i = length(angle_array):(-1):1
        % 注意这里的判断是 < d, 和之前的算法不同
        % 这里要求 range 的 min 角度对应距离 > d, max 角度对应距离 < d
        if func_dist(angel_1, angle_array(i)) > d   
            break
        end
    end

    % 牛顿迭代法获得高精度解

    % 进行牛顿迭代
    if i == 1
        min = angle_2_range(1);
        max = angle_array(1) + 0.1;
    else 
        min = angle_array(i);
        max = angle_array(i+1);
    end

    for i = 1:71     % ceil( (log(4*pi)+20*log(10))/log(2) ) = 71, 这使得 4*pi/2^n < 10^(-20)
        error_max = max - min;
        temp = 0.5 * (max + min);
        temp_dis = func_dist(angel_1, temp);
        if temp_dis < d
            max = temp;
        elseif temp_dis > d
            min = temp;
        elseif temp_dis == d
            angel = temp;
            return
        end
        if error_max <= error_max_abs
            break
        end
    end

    % 输出结果
    angel = 0.5*(min+max);
end

%% 问题一函数区 %%
%%%%%%%%%%%%%%%%%

function theta = GetNextPoint_newton( theta_0 , b, d , error_max_abs)
    % 第 1 层
    theta_range = theta_0:10^(-1):(theta_0 + 10*pi);
    for i = 2:length(theta_range)
        if GetDistance(theta_range(i), theta_0, b) > d
            break
        end
    end

    % 进行牛顿迭代
    min = theta_range(i-1);
    max = theta_range(i);

    for i = 1:71     % ceil( (log(4*pi)+20*log(10))/log(2) ) = 71, 这使得 4*pi/2^n < 10^(-20)
        error_max = max - min;
        temp = 0.5 * (max + min);
        temp_dis = GetDistance(temp, theta_0, b);
        if temp_dis < d
            min = temp;
        elseif temp_dis > d
            max = temp;
        elseif temp_dis == d
            theta = temp;
            return
        end
        if error_max <= error_max_abs
            break
        end
    end

    % 输出结果
    theta = 0.5*(min + max);
end


function dis = GetDistance(theta, theta_0, b) 
    dis = sqrt( b^2*(theta.^2 + theta_0.^2) -2*b^2.*theta.*theta_0.*cos(theta-theta_0) );
end


function theta_final = func_find_start(b, theta0, L)
    % 初始化参数
    % b = 55 / (2 * pi);    % 螺线递增参数 (b)
    % theta0 = 16 * 2 * pi; % 初始点的极坐标角度 (θ_0)
    % L = 500;              % 沿螺线走的距离 (L)
    
    % 定义螺旋线方程函数
    r = @(theta) b * (theta0-theta);       % 极径方程：r(θ) = b(θ_0-θ)
    ds = @(theta) sqrt(b^2 + r(theta).^2); % 曲线长度微元
    r_real = @(theta) b * theta;       % 极径方程：r(θ) = b(θ_0-θ)
    % 使用Numerical Integration计算θ值
    theta_final = fzero(@(theta) integral(ds, 0, theta) - L, theta0 + 2*pi);
    
end
