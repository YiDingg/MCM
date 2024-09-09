clear, clc, close all  
%% 程序主代码 %%
%%%%%%%%%%%%%%%%

% 数据准备
    d = 1.65;           % 板上节点距离 (m)
    d_0 = 2.86;         % 龙头板上节点距离 (m)
    luoju = 0.55;       % 螺距 (m)
    v_0 = 1;            % 龙头速度 (m/s)
    b = luoju/(2*pi);   % 计算 r = a + b*theta 的参数 b

% 等距螺线可视化
    figure
    stc.line = polarplot(0:0.1:(32*pi), b*(0:0.1:(32*pi)));
    hold on
    % 设置样式
        % 坐标轴
            stc.fig = gcf;
            stc.axes = gca;
            stc.axes.FontName = "Times New Roman"; % 全局 FontName
            stc.axes.Box = 'on';   
        % 图例
            stc.leg = legend(stc.axes, 'Location', 'northeast');
            stc.leg.FontSize = 15;
            stc.leg.Interpreter = "latex";
            stc.leg.String = '$r(\theta) = b\theta$';
            
        % 标题
            stc.axes.Title.String = '';
            stc.axes.Title.FontSize = 17;
            stc.axes.Title.FontWeight = 'bold';
        % 线的样式
            stc.line.LineWidth = 1;
            stc.line.Color = [0 0 1];   % 蓝色
        % 收尾
            hold(stc.axes,'off')
        %MyExport_pdf_docked

% 计算所有数据
    theta_array = zeros(224, 301);
    speed_array = zeros(224, 301);
    for t = 0:300
        theta_0 = 16*2*pi - func_find_start(b, 16*2*pi, v_0*t);
        theta_array(:, t+1) = GetAllPoints(theta_0, b, d, 10^(-16));
        speed_array(:, t+1) = GetAllSpeed(theta_array(:, t+1), b, v_0);
    end

% 输出结果
    PrintResult_Q1(theta_array, speed_array, b) 

% 部分结果可视化
    for t = [10, 290]
        theta_0 = 16*2*pi - func_find_start(b, 16*2*pi, v_0*t);
        theta_all = GetAllPoints(theta_0, b, d, 10^(-16));
        stc_drawpoints = DrawPoints(theta_all, b);
        stc_drawpoints.axes.Title.String = '';
        MyExport_pdf
        stc_rectangles = DrawPointsAndRectangles(theta_all, b);
        stc_rectangles.axes.Title.String = '';
        MyExport_pdf
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
    theta_range = linspace(theta_0, theta_0 + 10*pi, 300);
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
            stc.axes.FontSize = 11;
            xline(0, 'LineWidth', 0.3, 'Color',  [0.7, 0.7, 0.7]);
            yline(0, 'LineWidth', 0.3, 'Color',  [0.7, 0.7, 0.7]);
            stc.label.x = xlabel(stc.axes, '$x\ (\mathrm{m})$', 'Interpreter', 'latex', 'FontSize', 15);
            stc.label.y = ylabel(stc.axes, '$y\ (\mathrm{m})$', 'Interpreter', 'latex', 'FontSize', 15);
            

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
        %scatter(Rectangle_Points(1, i, 1), Rectangle_Points(1, i, 2), 45, '.red');
        %scatter(Rectangle_Points(2, i, 1), Rectangle_Points(2, i, 2), 45, '.black');
        %scatter(Rectangle_Points(3, i, 1), Rectangle_Points(3, i, 2), 45, '.black');
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
