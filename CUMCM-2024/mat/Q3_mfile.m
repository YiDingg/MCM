clc, clear, close all
%% 主程序代码区 %%
%%%%%%%%%%%%%%%%%

% 模拟退火获取最小螺距
    
    % 数据准备
        d = 1.65;           % 板上节点距离 (m)
        d_0 = 2.86;         % 龙头板上节点距离 (m)
        v_0 = 1;            % 龙头速度 (m/s)
        t_range = [0, 350];
        t_step1 = 0.5;
        t_step2 = 0.1;
        t_error_max_abs = 10^(-12);
    
    % 构建退火结构体
        Struct_SA.Var.num = 1; % 一个参数
        Struct_SA.Var.range = [0.44 0.46 0.455];    % 螺距参数范围和初始值
        Struct_SA.Annealing.T_0 = 10;        % 初始温度
        Struct_SA.Annealing.alpha = 0.97;    % 降温系数
        Struct_SA.Annealing.T_end = 1;       % 停止温度
        Struct_SA.Annealing.mkvlength = 4;   % 马尔科夫链长度

% 进行模拟退火
    Struct_SA = MySimulatedAnnealing(Struct_SA, @Objective_Annealing, d, v_0, t_range, t_step1, t_step2, t_error_max_abs);
    vpa(Struct_SA.X_best)   % 显示结果


% 可靠性检验
    % 图 1
    luoju_range = [0.35 0.5];
    luoju_array = linspace(luoju_range(1), luoju_range(2), 50);
    dist_array = zeros(size(luoju_array));
    for i = 1:length(luoju_array)
        dist_array(i) = GetCurushedDistance(luoju_array(i), d, v_0, t_range, t_step1, t_step2, t_error_max_abs);
    end
    stc = MyPlot(luoju_array, dist_array);
    yline(4.5);
    stc.label.x.String = 'screw pitch (m)';
    stc.label.y.String = 'distance (m)';
    stc.axes.Title.String = '';
    stc.leg.String = ['crushing distance'; ''; ''; ''; ''; ''];

    % 图 2
    luoju_range = [0.44 0.46];
    luoju_array = linspace(luoju_range(1), luoju_range(2), 50);
    dist_array = zeros(size(luoju_array));
    for i = 1:length(luoju_array)
        dist_array(i) = GetCurushedDistance(luoju_array(i), d, v_0, t_range, t_step1, t_step2, t_error_max_abs);
    end
    stc = MyPlot(luoju_array, dist_array);
    yline(4.5);
    stc.label.x.String = 'screw pitch (m)';
    stc.label.y.String = 'distance (m)';
    stc.axes.Title.String = '';
    stc.leg.String = ['crushing distance'; ''; ''; ''; ''; ''];

%% 问题三函数区 %%
%%%%%%%%%%%%%%%%%

function obj = Objective_Annealing(luoju, d, v_0, t_range, t_step1, t_step2, t_error_max_abs)
    dist = GetCurushedDistance(luoju, d, v_0, t_range, t_step1, t_step2, t_error_max_abs);
    if dist < 4.5   % 能进入掉头区
        obj = luoju;
    else    % 不能进入调头区
        obj = 0.455;
    end
end 

function dist = GetCurushedDistance(luoju, d, v_0, t_range, t_step1, t_step2, t_error_max_abs)
% 函数：获取最小螺距
    b = luoju/(2*pi);   % 计算 r = a + b*theta 的参数 b
    t_crushed = GetCrushedTime(b, d, v_0, t_range, t_step1, t_step2, t_error_max_abs);
    theta_0 = 16*2*pi - func_find_start(b, 16*2*pi, v_0*t_crushed);
    coordinates = b*theta_0.*[cos(theta_0), sin(theta_0)];
    dist = sqrt(sum(coordinates.^2));
    %disp(['luoju = ', num2str(luoju), ', t_crushed = ', num2str(t_crushed), ', theta_0 = ', num2str(theta_0) ', distance = ', num2str(dist) ])
    %DrawPointsAndRectangles(GetAllPoints(theta_0 , b, d , 10^(-12)), b);
end

%% 问题二函数区 %%
%%%%%%%%%%%%%%%%%

function t_crushed = GetCrushedTime(b, d, v_0, t_range, t_step1, t_step2, error_max_abs)
    % 第一层粗搜
    for t = t_range(1):t_step1:t_range(2)
        theta_0 = 16*2*pi - func_find_start(b, 16*2*pi, v_0*t);
        theta_all = GetAllPoints(theta_0, b, d, 10^(-16));
        logic = IsCrushed(theta_all, b);
        if logic == true
            break
        end
    end

    % 第二层细搜
    t_range = [t - 3*t_step1, t];
    for t = t_range(1):t_step2:t_range(2)
        theta_0 = 16*2*pi - func_find_start(b, 16*2*pi, v_0*t);
        theta_all = GetAllPoints(theta_0, b, d, 10^(-16));
        logic = IsCrushed(theta_all, b);
        if logic == true
            break
        end
    end

    % 牛顿迭代法获取高精度解
    min = t - t_step2;  max = t;
    for i = 1:67     % ceil( (log(1)+20*log(10))/log(2) ) = 67, 这使得 1/2^n < 10^(-20)
        error_max = max - min;
        temp_t = 0.5 * (max + min);
        theta_0 = 16*2*pi - func_find_start(b, 16*2*pi, v_0*temp_t);
        theta_all = GetAllPoints(theta_0, b, d, 10^(-16));
        logic = IsCrushed(theta_all, b);
        if logic ~= true
            min = temp_t;
        elseif logic == true
            max = temp_t;
        end
        if error_max <= error_max_abs
            break
        end
    end

    % 输出结果
    t_crushed = 0.5*(min+max);
end

function logic = IsCrushed(theta_all, b)
    % 数据准备
    d_special = sqrt(3.135^2 + 0.15^2) + sqrt(0.275^2 + 0.15^2);  % 龙头 3 号点最大判断间距
    d_common = sqrt(1.925^2 + 0.15^2) + sqrt(0.275^2 + 0.15^2);  % 2, 3 号点最大判断间距
    Rectangle_Points = GetRectangles(theta_all, b);

    % 只需判断第 2, 3 号点是否 legal

    % 判断龙头 (第一个板凳)
        point_1_2(1, 1:2) = Rectangle_Points(2, 1, :);
        point_1_3(1, 1:2) = Rectangle_Points(3, 1, :);

        % 挑选可能的节点
        distance = sqrt( b^2*(theta_all(1).^2 + theta_all.^2) -2*b^2.*theta_all(1).*theta_all.*cos(theta_all(1)-theta_all) );
        Logic = (distance < d_special) ;
        Logic(1:2) = 0; % 忽略第 1(本身),2(之后)
        Logic(224) = 0; % 忽略最后一个节点 (无板凳)
        for i = find(Logic)'    % 这里必须转置使得右值为行向量
            % 新的坐标原点
            new_Origin(1, 1:2) = Rectangle_Points(1, i, :);

            % 检查 2 号点 
            new_cor = (point_1_2 - new_Origin ) / [
                Rectangle_Points(2, i, 1) - Rectangle_Points(1, i, 1), Rectangle_Points(2, i, 2) - Rectangle_Points(1, i, 2);
                Rectangle_Points(4, i, 1) - Rectangle_Points(1, i, 1), Rectangle_Points(4, i, 2) - Rectangle_Points(1, i, 2);
                ];   % 坐标系转换
            %disp(['龙头2号点 new_cor:', num2str(new_cor)])
            if (0<new_cor(1)&&new_cor(1)<1) && (0<new_cor(2)&&new_cor(2)<1)
                logic = true;
                return
            end
            % 检查 3 号点
            new_cor = (point_1_3 - new_Origin ) / [
                Rectangle_Points(2, i, 1) - Rectangle_Points(1, i, 1), Rectangle_Points(2, i, 2) - Rectangle_Points(1, i, 2);
                Rectangle_Points(4, i, 1) - Rectangle_Points(1, i, 1), Rectangle_Points(4, i, 2) - Rectangle_Points(1, i, 2);
                ];   % 坐标系转换
            %disp(['龙头3号点 new_cor:', num2str(new_cor)])
            if (0<new_cor(1)&&new_cor(1)<1) && (0<new_cor(2)&&new_cor(2)<1)
                logic = true;
                return
            end
        end


    % 判断第二个板凳
        point_2_2(1, 1:2) = Rectangle_Points(2, 2, :);
        point_2_3(1, 1:2) = Rectangle_Points(3, 2, :);

        % 挑选可能的节点
        distance = sqrt( b^2*(theta_all(2).^2 + theta_all.^2) -2*b^2.*theta_all(2).*theta_all.*cos(theta_all(2)-theta_all) );
        Logic = (distance < d_common) ;
        Logic(1:3) = 0; % 忽略第 1(龙头),2(本身),3(其后) 节点
        Logic(224) = 0; % 忽略最后一个节点 (无板凳)
        for i = find(Logic)'    % 这里必须转置使得右值为行向量
            % 新的坐标原点
            new_Origin(1, 1:2) = Rectangle_Points(1, i, :);

            % 检查 2 号点
            new_cor = (point_2_2 - new_Origin ) / [
                Rectangle_Points(2, i, 1) - Rectangle_Points(1, i, 1), Rectangle_Points(2, i, 2) - Rectangle_Points(1, i, 2);
                Rectangle_Points(4, i, 1) - Rectangle_Points(1, i, 1), Rectangle_Points(4, i, 2) - Rectangle_Points(1, i, 2);
                ];   % 坐标系转换
            %disp(['板凳2号点 new_cor:', num2str(new_cor)])
            if (0<new_cor(1)&&new_cor(1)<1) && (0<new_cor(2)&&new_cor(2)<1)
                logic = true;
                return
            end
            % 检查 3 号点
            new_cor = (point_2_3 - new_Origin) / [
                Rectangle_Points(2, i, 1) - Rectangle_Points(1, i, 1), Rectangle_Points(2, i, 2) - Rectangle_Points(1, i, 2);
                Rectangle_Points(4, i, 1) - Rectangle_Points(1, i, 1), Rectangle_Points(4, i, 2) - Rectangle_Points(1, i, 2);
                ];   % 坐标系转换
            %disp(['板凳3号点 new_cor:', num2str(new_cor)])
            if (0<new_cor(1)&&new_cor(1)<1) && (0<new_cor(2)&&new_cor(2)<1)
                logic = true;
                return
            end
        end

    % 检查通过
    logic = false;
    return
end


function PrintResult_Q2(theta, Speed, b)   
    % theta: 224*1 矩阵
    % Speed: 224*1 矩阵
    % b: 用于转直角坐标
    % 注：未进圈的数据本应是 nan, 为方便考虑 
    %     这里未做处理, theta_all((theta_all - 16*2*pi) > 0) = nan 即可筛选

    % 计算直角坐标位置
    Coordinates = b*theta.*[cos(theta), sin(theta)];

    % 遍历矩阵，将每个元素格式化为保留 6 位小数的字符串  
    Output = num2cell(zeros(224,3));
    for i = 1:size(Coordinates, 1)  
            Output{i, 1} = num2str(Coordinates(i, 1), '%.6f');  
            Output{i, 2} = num2str(Coordinates(i, 2), '%.6f');  
    end  
    for i = 1:size(Speed, 1)  
        Output{i, 3} = num2str(Speed(i), '%.6f');  
    end  

    writecell(Output, 'Q2_Result.xlsx', 'Sheet', 'Sheet1'); % 输出结果
end


%% 问题一函数区 %%
%%%%%%%%%%%%%%%%%

function dis = GetDistance(theta, theta_0, b) 
    dis = sqrt( b^2*(theta.^2 + theta_0.^2) -2*b^2.*theta.*theta_0.*cos(theta-theta_0) );
end

% 此函数已废弃，之后使用牛顿迭代法
function [theta, rho] = GetNextPoint( theta_0 , b, d , error_level)
    % 第 1 层
    theta_range = theta_0:10^(-1):(theta_0 + 10*pi);
    for i = 2:length(theta_range)
        if GetDistance(theta_range(i), theta_0, b) > d
            break
        end
    end

    % 第 2 ~ error_level 层
    for j = 2:error_level
        theta_range = theta_range(i-1):10^(-j):theta_range(i);
        for i = 2:length(theta_range)
            if GetDistance(theta_range(i), theta_0, b) > d
                break
            end
        end
    end

    % 输出结果
    theta = 0.5*( theta_range(i-1) + theta_range(i) );
    rho = b*theta;
end

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

function theta_all = GetAllPoints( theta_0 , b, d , error_max_abs)
    theta_all = zeros(224, 1);
    theta_all(1) = theta_0;
    % 龙头板凳长不同
    theta_all(2) = GetNextPoint_newton(theta_all(1), b, 2.86, error_max_abs);
    for i = 3:224
        theta_all(i) = GetNextPoint_newton(theta_all(i-1), b, d, error_max_abs);
        %polarscatter(theta, r, 100, 'r.');
    end
end


function stc = DrawPoints(theta_all, b)

    X = 0:0.1:(32*pi);
    Rho = b*X;

    % 作图
    figure
    stc.line = polarplot(X, Rho);
    hold on
    stc.point_0 = polarscatter(theta_all(1), b*theta_all(1), 200, 'r.');
    stc.points_rest = polarscatter(theta_all(2:end), b*theta_all(2:end), 70, 'b.');

% 设置样式
    % 坐标轴
        stc.fig = gcf;
        stc.axes = gca;
        stc.axes.FontName = "Times New Roman"; % 全局 FontName
        stc.axes.Box = 'on';  

    % 标题
        stc.axes.Title.String = 'Figure: Draw Points';
        stc.axes.Title.FontSize = 17;
        stc.axes.Title.FontWeight = 'bold';
    
    % 线的样式
        stc.line.LineWidth = 1;
        stc.line.Color = [0 1 0];   % 绿色

    % 收尾
        hold(stc.axes,'off')
end


function Rectangle_Points = GetRectangles(theta_all, b)
    Coordinates = b*theta_all.*[cos(theta_all), sin(theta_all)];
    Vec_X = diff(Coordinates);
    Vec_X = Vec_X ./ sqrt(sum(Vec_X.^2, 2));     % 单位化
    Vec_N = [ -Vec_X(: ,2), Vec_X(: ,1) ];  % 法向量


    % 计算矩形坐标
    Rectangle_P1 = Coordinates(1:223, :) - Vec_X*0.275 + Vec_N*0.15;    % 注意是 1:223
    Rectangle_P2 = Coordinates(1:223, :) - Vec_X*0.275 - Vec_N*0.15;
    Rectangle_P3 = Coordinates(2:224, :) + Vec_X*0.275 - Vec_N*0.15;    % 注意是 2:224
    Rectangle_P4 = Coordinates(2:224, :) + Vec_X*0.275 + Vec_N*0.15;

    Rectangle_Points = zeros(4, 223, 2);
    Rectangle_Points(1, :, :) = Rectangle_P1;
    Rectangle_Points(2, :, :) = Rectangle_P2;
    Rectangle_Points(3, :, :) = Rectangle_P3;
    Rectangle_Points(4, :, :) = Rectangle_P4;
end


function stc = DrawPointsAndRectangles(theta_all, b)
% 作线和点
    x = 0:0.1:(32*pi);
    X = x';

    % 转为直角坐标
        coor_line = [
            b*X.*cos(X), b*X.*sin(X)
        ];
        coor_all = [
            b*theta_all.*cos(theta_all), b*theta_all.*sin(theta_all)
        ];

    % 作图
        figure
        stc.line = plot(coor_line(:, 1), coor_line(:, 2));
        hold on
        stc.point_0 = scatter(coor_all(1, 1), coor_all(1, 2), 200, 'r.');
        stc.points_rest = scatter(coor_all(2:end, 1), coor_all(2:end, 2), 70, 'b.');

    % 设置样式
        % 坐标轴
            stc.fig = gcf;
            axis equal  
            stc.axes = gca;
            stc.axes.FontName = "Times New Roman"; % 全局 FontName
            stc.axes.Box = 'on';  
            stc.axes.FontSize = 14;
            xline(0, 'LineWidth', 0.3, 'Color',  [0.3, 0.3, 0.3]);
            yline(0, 'LineWidth', 0.3, 'Color',  [0.3, 0.3, 0.3]);
            stc.label.x = xlabel(stc.axes, '$x$', 'Interpreter', 'latex', 'FontSize', 15);
            stc.label.y = ylabel(stc.axes, '$y$', 'Interpreter', 'latex', 'FontSize', 15);
            

        % 标题
            stc.axes.Title.String = 'Figure: Draw Rectangles';
            stc.axes.Title.FontSize = 17;
            stc.axes.Title.FontWeight = 'bold';
        
        % 线的样式
            stc.line.LineWidth = 1;
            stc.line.Color = [0 1 0];   % 绿色

% 作方框 
    Rectangle_Points = GetRectangles(theta_all, b);
    Rectangle_Points(5, :, :) = Rectangle_Points(1, :, :);    % plot 围成闭合曲线
    hold on
    for i = 1:223
        plot(Rectangle_Points(:, i, 1), Rectangle_Points(:, i, 2), 'LineWidth', 0.3, 'Color', [1 0 1]);
    end

% 收尾
    hold off
end


function Speed_all = GetAllSpeed(theta_all, b, v_0)
    Speed_all = zeros(224 ,1);
    Speed_all(1) = v_0;

    % 转为直角坐标
    Coordinates = b*theta_all.*[cos(theta_all), sin(theta_all)];
    Vec_X = diff(Coordinates);
    % 下面注意要有负号
    Vec_Tao = - [cos(theta_all) - theta_all.*sin(theta_all),  sin(theta_all) + theta_all.*cos(theta_all)];
    for i = 1:223
        Speed_all(i+1) = Speed_all(i) * ( norm(Vec_Tao(i+1, :))*abs(sum(Vec_Tao(i, :).*Vec_X(i, :)))  ) / ( norm(Vec_Tao(i, :))*abs(sum(Vec_Tao(i+1, :).*Vec_X(i, :))) );
    end
end


% 此函数已废弃，因为 .xlsx 格式不同
function printLocation(theta_all, Speed_all, b, column_initial)
    %logic_matrix = (theta_all - 16*2*pi) > 0;
    theta_all((theta_all - 16*2*pi) > 0) = nan;    % 将不合法 (未进圈) 的数据设为 nan 
    Speed_all((theta_all - 16*2*pi) > 0) = nan;
    % 计算直角坐标位置
    Coordinates = b*theta_all.*[cos(theta_all), sin(theta_all)];
    rowHeaders = arrayfun(@(i) sprintf('第%d节点', i), 1:length(theta_all), 'UniformOutput', false); 

    Output = array2table([Coordinates, Speed_all], "RowNames",rowHeaders, "VariableNames", {'位置 x', '位置 y', '速率 v'});
    writetable(Output, 'Q1_Result.xlsx', 'WriteRowNames', true, 'Range', column_initial);
end

function PrintResult_Q1(theta, Speed, b)   
    % theta: 224*301 矩阵
    % Speed: 224*301 矩阵
    % b: 用于转直角坐标
    % 注：未进圈的数据本应是 nan, 为方便考虑 
    %     这里未做处理, theta_all((theta_all - 16*2*pi) > 0) = nan 即可筛选

    t_length = size(theta, 2);

    % 计算直角坐标位置
    Coordinates = zeros(2, 224, t_length);
    Coordinates(1, :, :) = b*theta.*cos(theta);
    Coordinates(2, :, :) = b*theta.*sin(theta);
    Location = zeros(2*224, t_length);
    for i = 1:2*224
        if mod(i, 2) == 1   % i 奇数, 对应 x 
            Location(i, :) = Coordinates(1, ceil(i/2), :);
        else                % i 偶数, 对应 y 
            Location(i, :) = Coordinates(2, i/2, :);
        end
    end
    % 遍历矩阵，将每个元素格式化为保留 6 位小数的字符串  
    for i = 1:size(Location, 1)  
        for j = 1:size(Location, 2)  
            Location_str{i, j} = num2str(Location(i, j), '%.6f');  
        end  
    end  
    for i = 1:size(Speed, 1)  
        for j = 1:size(Speed, 2)  
            Speed_str{i, j} = num2str(Speed(i, j), '%.6f');  
        end  
    end  

    writecell(Location_str, 'Q1_Result.xlsx', 'Sheet', 'Sheet1'); % 输出位置
    writecell(Speed_str, 'Q1_Result.xlsx', 'Sheet', 'Sheet2'); % 输出速度
end

function theta = GetTheta(b, v, t)
    % 弧长公式
    S = @(theta) b/2 * (  theta.*sqrt(1+theta.^2) + log(theta + sqrt(1+theta.^2))  );
    S_real = @(theta) S(16*2*pi) - S(theta);
    
    % 第 1 层
    theta_range = 16*2*pi:(-10^(-1)):5*pi;
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
    % 输入参数
    % b = 55 / (2 * pi);    % 螺线递增参数 (b)
    % theta0 = 16 * 2 * pi; % 初始点的极坐标角度 (θ_0)
    % L = 500;              % 沿螺线走的距离 (L)
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

