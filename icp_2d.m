%% SVD 解析解的 2D ICP 算法实现
function icp_2d()
    visuInit();

    name = "laserscan.csv";
    lasercan = csvread(name);

    IcpAndVius(lasercan, 4.5, 0.33);
    IcpAndVius(lasercan, 4.5, 12);
    IcpAndVius(lasercan, 16, 0.33);
    IcpAndVius(lasercan, 16, 12);

end

function IcpAndVius(scan, init_angle, init_offset)
    % 对 scan 做 init_angle, init_offset 给定的偏移，作为待配准的点云
    theta = init_angle * pi / 180;
    R = [
        cos(theta), -sin(theta);
        sin(theta), cos(theta);
    ];
    t = [init_offset; 0];
    
    X = scan';
    Y = R * X + t;

    % ICP 迭代
    for i = 1:16
        visu(X, Y, init_angle, init_offset, i);
        if i == 1
            pause(3);
        else 
            pause(1);
        end

        [R, t] = computeTransform(X, Y);
        Y = R * Y + t;
    end
end

% 求点云的质心
function result = centroid(point_cloud)
    point_count = size(point_cloud, 2);
    result = sum(point_cloud, 2) / point_count;
end

% SVD 解旋转矩阵和位移
% X 目标点云，规模 3xN 的 2d 齐次坐标
% Y 待配准点云，规模 3xN 的 2d 齐次坐标
% [R, t] Y 对其 X 需要的线性变换
function [R, t] = computeTransform(X, Y)
    knn_index = knnsearch(X', Y');
    disp("knn_index");
    disp(size(knn_index, 1))

    % 质心
    X_centroid = centroid(X);
    Y_centroid = centroid(Y);
    
    X_residual = X - X_centroid;
    Y_residual = Y - Y_centroid;

    % H = Y_residual * X_residual';
    % 临近点构造 H
    H = [
        0 0; 
        0 0;
    ];
    for i = 1:size(knn_index, 1)
        index_y = knn_index(i);
        H = H + (Y_residual(:, index_y) * X_residual(:, i)');
    end

    [U, S, V] = svd(H);
    R = V * U';
    t = X_centroid - R * Y_centroid;
end



%% 可视化
function visuInit()
    fig = figure;
    pos = fig.Position; 
    pos(1) = 0;
    pos(2) = 0;
    pos(3) = 800; 
    pos(4) = 600; 
    fig.Position = pos; 
end

function visu(X, Y, init_angle, init_offset, i)
    X_centroid = centroid(X);
    Y_centroid = centroid(Y);
    hold off;
    scatter(X(1, :), X(2, :), 8, 'filled', 'MarkerFaceColor', [0, 0, 1], 'DisplayName', 'X');
    hold on;
    scatter(X_centroid(1, :), X_centroid(2, :), 30, 'filled', 'MarkerFaceColor', [0, 0.2, 1], 'DisplayName', 'X centroid');
    
    scatter(Y(1, :), Y(2, :), 4, 'filled', 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'Y');
    scatter(Y_centroid(1, :), Y_centroid(2, :), 20, 'filled', 'MarkerFaceColor', [1, 0.2, 0], 'DisplayName', 'Y centroid');

    axis equal;
    xlabel("x");
    ylabel("y");
    title(sprintf("初始角度差 %.3f °，初始距离差 %.3f m\n迭代 %d 次",init_angle, init_offset, i));
    grid on;
    daspect([1 1 1]);
    xlim([-10, 25]);
    ylim([-10, 20]);
end
