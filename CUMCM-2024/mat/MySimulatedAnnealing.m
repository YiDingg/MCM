function [stc_SA, stc_Figure] = MySimulatedAnnealing(stc_SA, objective, d, v_0, t_range, t_step1, t_step2, t_error_max_abs)
% 输入退火问题结构体，输出迭代结果（minimize）。
% 
% 输入：
    % stc_SA：退火问题结构体
    % objective：目标函数
% 输出：迭代结果
% 注：迭代总次数为 1000 时，在 waitbar 上共耗时约 0.8 s 
%%
    % 步骤一：初始化
        %Waitbar = waitbar(0, '1', 'Name', 'Simulated Annealing', 'CreateCancelBtn', 'delete(gcbf)','Color', [0.9, 0.9, 0.9]);
        %Btn = findall(Waitbar, 'type', 'Uicontrol');
        %set(Btn, 'String', 'Cancle', 'FontSize', 10);
        Waitbar = waitbar(0, '1', 'Name', 'Simulated Annealing', 'Color', [0.9, 0.9, 0.9]);
        TK = stc_SA.Annealing.T_0;
        T_end = stc_SA.Annealing.T_end;
        mkv = stc_SA.Annealing.mkvlength;
        alpha = stc_SA.Annealing.alpha;
        num_Var = size(stc_SA.Var.range(:,3));
        num_Var = num_Var(1);
        %lam = stc.Annealing.lam;
        % 迭代记录仪初始化
        change_1 = 0;    % 随机到更优解的次数
        change_2 = 0;    % 接收较差解的次数（两者约 10:1 时有较好寻优效果）
        mytry = 1;       % 当前迭代次数
        N = mkv*ceil(log(T_end/TK)/log(alpha));   % 迭代总次数
        X = stc_SA.Var.range(:,3)';
        X_best = stc_SA.Var.range(:,3)';
        f_best = objective(X, d, v_0, t_range, t_step1, t_step2, t_error_max_abs); % 计算目标函数
        process = zeros(1, floor(N));
        process(1) = f_best;
        process_best = zeros(1, floor(N));
        process_best(1) = f_best;
        

    % 步骤二：退火
        disp("初始化完成，开始退火")
        disp(['initial obj: ',num2str(f_best)])
        start = tic; % 开始计时
        while TK >= T_end   
            for i = 1:mkv  % 每个温度T下，我们都寻找 mkv 次新解 X，每一个新解都有可能被接受

                r = rand;
                if r>=0.5 % 在当前较优解附近扰动
                    for j = 1:num_Var
                            X(j) = X_best(j)+(rand-0.5)*(stc_SA.Var.range(j,2) - stc_SA.Var.range(j,1))*( 1-(mytry-1)/N )^2;
                            X(j) = max(stc_SA.Var.range(j,1), min(X(j), stc_SA.Var.range(j,2)));   % 确保扰动后的 X 仍在范围内
                    end
                else % 生成全局随机解
                    X = (stc_SA.Var.range(:,1) + rand*(stc_SA.Var.range(:,2) - stc_SA.Var.range(:,1)))';  % 转置后才是行向量   
                end

                mytry = mytry+1;
                f = objective(X, d, v_0, t_range, t_step1, t_step2, t_error_max_abs); % 计算目标函数

                if f < f_best   % 随机到更优解，接受新解
                   f_best = f;
                   X_best = X;
                   change_1 = change_1+1;
                   % disp(['较优参数为：',num2str(X_best)])
                   disp(['    new obj: ',num2str(f_best)])
                elseif exp( - 10^7 * abs((f_best-f)/f_best) /TK  ) > rand  % 满足概率，接受较差解
                   f_best = f;
                   X_best = X;
                   % disp(['较优参数为：',num2str(X_best)])
                   disp(['    new obj: ',num2str(f_best)])
                   change_2 = change_2 + 1;
                end

                process(mytry) = f;
                process_best(mytry) = f_best;
            end
            %disp(['进度：',num2str((mytry)/N*100),'%'])
            TK = TK*alpha;
            waitbar((mytry)/N*100, Waitbar, ['Computing: ', num2str(round((mytry)/N*100, 1)), '%']);
        end

    % 步骤三：退火结束，输出最终结果
        % 进度条
            time = toc(start);
            waitbar(1, Waitbar, ['Simulated Annealing Completed (in ', num2str(time),' s)']);
            Waitbar.Color = [1 1 1];
        % 图像
            stc_SA.process = process;
            stc_SA.process_best = process_best;
            stc_SA.X_best = X_best;
            stc_SA.Object_best = f_best;

            
            stc_Figure = MyYYPlot(1:length(process), process_best, 1:length(process), process);
            stc_Figure.axes.Title.String = ['Simulated Annealing (in ', num2str(time),' s)'];
            %stc_Figure.axes.YLimitMethod = "padded";
            stc_Figure.leg.String = ["Best objective value";"Current objective value"];
            stc_Figure.leg.Location = "northeast";
            stc_Figure.p_left.LineWidth = 5;
            stc_Figure.p_right.LineWidth = 1;
            stc_Figure.label.x.String = 'times';
            stc_Figure.label.y_left.String = '$obj_{\mathrm{best}}$';
            stc_Figure.label.y_right.String = '$obj_{\mathrm{current}}$';

        % 文本
            disp('---------------------------------')
            disp('>> --------  模拟退火  -------- <<')
            disp(['历时 ', num2str(time), ' 秒'])
            disp(['一共寻找新解：',num2str(mytry)])
            disp(['change_1次数：',num2str(change_1)])
            disp(['change_2次数：',num2str(change_2)])
            disp(['最优参数：', num2str(X_best)])
            disp(['最优目标值：', num2str(f_best)])
            disp('>> --------  模拟退火  -------- <<')
            disp('---------------------------------')

        % 导出数据
            output = cell(1, 6 + num_Var);
            output{1, 1} = datestr(now, 'yyyy-mm-dd HH:MM:SS');
            output{1, 2} = 'Spend (s)';
            output{1, 3} = num2str(time, '%.6f');
            output{1, 4} = 'X_best';
            for i = 1:num_Var
                output{1, 4+i} = num2str(stc_SA.X_best(i), '%.12f');
            end
            output{1, 5+num_Var} = 'Object_best';
            output{1, 6+num_Var} = num2str(stc_SA.Object_best, '%.12f');
            writecell(output, 'MySimulatedAnnealingResualts.xlsx', "WriteMode","append");
            %writematrix([datestr(now, 'yyyy-mm-dd HH:MM:SS'), "Spend (s)", time, 'X_best', vpa(stc_SA.X_best), 'Object_best', vpa(stc_SA.Object_best)], 'MySimulatedAnnealingResualts.xlsx', "WriteMode","append");
end
