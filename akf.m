function [x, P] = akf(x,P,F,H,z,Q,R,lambda)
    % 自适应卡尔曼滤波器
    % z: 测量值向量
    % Q: 过程噪声协方差
    % R: 初始测量噪声协方差
    % lambda: 自适应调整因子 (0 < lambda < 1)
    
    % 时间更新 (预测)
    x_pred = F*[x(1); x(1:n-1)]; 
    P_pred = F*P*F' + G*Q*G'; 

    % 测量更新 (校正)
    y = z(k) - H*x_pred; % 测量残差
    S = P_pred + R; % 测量残差协方差
    K = P_pred / S; % 卡尔曼增益
    x(k) = x_pred + K * y; % 状态更新
    P(k) = (1 - K) * P_pred; % 协方差更新

    % 自适应调整测量噪声协方差
    R = (1 - lambda) * R + lambda * S; % 指数平滑调整 R

end