function attitude = ins_attitude(gyro,acc,init_yaw)

g = 9.7970;

% %陀螺零偏
% gyrox_bais = mean(gyro(1:1500,1));
% gyroy_bais = mean(gyro(1:1500,2));
% gyroz_bais = 11.5*mean(gyro(1:1000,3));
% % gyroz_bais = 11.5*mean(gyro(1:1500,3));
% 
% %加计零偏
% accx_bais = mean(acc(1:1000,1));
% accy_bais = mean(acc(1:1000,2));
% accz_bais = mean(acc(1:1000,3)-g);

%陀螺零偏
gyrox_bais = mean(gyro(1:50,1));
gyroy_bais = mean(gyro(1:50,2));
gyroz_bais = 11.5*mean(gyro(1:50,3));
% gyroz_bais = 11.5*mean(gyro(1:1500,3));

%加计零偏
accx_bais = mean(acc(1:50,1));
accy_bais = mean(acc(1:50,2));
accz_bais = mean(acc(1:50,3)-g);

%参数设置
n_start = 2;%实测数据开始数
leng = length(acc(:,1));%数组长度

% ACC_K_ = [2.768486046581383e-06,-1.289944236627479e-08,2.336544303499619e-08,
%           2.165439592511415e-08,2.773246582199320e-06,-4.393657219168412e-08,
%           2.693251278847677e-08,3.112040419592376e-08,2.755631130574960e-06];
% ACC_K0 = [4.209233997664602e+06,4.258163073880744e+06,4.195638544710741e+06];
% Kg_ = [1.00297946863758,1.00103643070260,1.00520913930355];
% Kg0 = [0.129020652898068,-0.036994003997335,-0.096067754830113];

% g = 9.7997;

    fs = 100;

%     Kp = 0.5;    
%     Ki = 0.005;                                        
% Kp = 0;    
% Ki = 0; 
      
    halfT = 1/fs/2;
    exInt(n_start) = 0;
    eyInt(n_start) = 0;
    ezInt(n_start) = 0;
    
    norm(1) = sqrt(acc(1,1)*acc(1,1) + acc(1,2)*acc(1,2) + acc(1,3)*acc(1,3));

    attitude(1,1) = asin(acc(1,2)/norm(1));%pitch俯仰 
    attitude(1,2) = -atan2(acc(1,1)/norm(1),acc(1,3)/norm(1));%roll翻滚
    attitude(1,3) = (init_yaw-180)*0.0174533;%此处对准完成   atan2(y, x)

    si = sin(attitude(1,1)); ci = cos(attitude(1,1)); %俯仰pit
    sj = sin(attitude(1,2)); cj = cos(attitude(1,2)); %翻滚roll
    sk = sin(attitude(1,3)); ck = cos(attitude(1,3)); %航向head
    
    Cnb = [       cj*ck,       cj*sk,        -sj;
            si*sj*ck-ci*sk, si*sj*sk+ci*ck, si*cj;
            ci*sj*ck+si*sk, ci*sj*sk-si*ck, ci*cj ];
   %P180_7姿态阵转化为四元数     
    q0(n_start) = sqrt(abs(1.0 + Cnb(1,1) + Cnb(2,2) + Cnb(3,3)))/2.0;
    q1(n_start) = sign(Cnb(3,2)-Cnb(2,3)) * sqrt(abs(1.0 + Cnb(1,1) - Cnb(2,2) - Cnb(3,3)))/2.0;
    q2(n_start) = sign(Cnb(1,3)-Cnb(3,1)) * sqrt(abs(1.0 - Cnb(1,1) + Cnb(2,2) - Cnb(3,3)))/2.0;
    q3(n_start) = sign(Cnb(2,1)-Cnb(1,2)) * sqrt(abs(1.0 - Cnb(1,1) - Cnb(2,2) + Cnb(3,3)))/2.0; %sign判断正负，<0为-1，>0为1，=0为0
    
    normq(n_start-1) = sqrt(q0(n_start)*q0(n_start) + q1(n_start)*q1(n_start) + q2(n_start)*q2(n_start) + q3(n_start)*q3(n_start));
    q0(n_start) = q0(n_start) / normq(n_start-1);
    q1(n_start) = q1(n_start) / normq(n_start-1);
    q2(n_start) = q2(n_start) / normq(n_start-1);
    q3(n_start) = q3(n_start) / normq(n_start-1);
for i=n_start:leng

    gx(i) = (gyro(i,1)-gyrox_bais) * 0.0174533; %0.0174533角度转弧度
    gy(i) = (gyro(i,2)-gyroy_bais) * 0.0174533;
    gz(i) = (gyro(i,3)-gyroz_bais) * 0.0174533;

    acx(i) = acc(i,1)-accx_bais;
    acy(i) = acc(i,2)-accy_bais;
    acz(i) = acc(i,3)-accz_bais;
    
    if(~((acx(i) == 0) && (acy(i) == 0) && (acz(i) == 0))) 
        %第一步重力加速度归一化
        norm(i) = sqrt(acx(i)*acx(i) + acy(i)*acy(i) + acz(i)*acz(i));
        ax(i) = acx(i) / norm(i);
        ay(i) = acy(i) / norm(i);
        az(i) = acz(i) / norm(i);
        
        %第二步提取四元数的等效余弦矩阵中的重力分量，陀螺
        vx(i) = 2*(q1(i)*q3(i) - q0(i)*q2(i));	
        vy(i) = 2*(q0(i)*q1(i) + q2(i)*q3(i));
        vz(i) = q0(i)*q0(i) - q1(i)*q1(i) - q2(i)*q2(i) + q3(i)*q3(i);		
	        
        %第三步向量叉积得出姿态误差    把矩阵写出来
        ex(i) = (ay(i)*vz(i) - az(i)*vy(i));
        ey(i) = (az(i)*vx(i) - ax(i)*vz(i));
        ez(i) = (ax(i)*vy(i) - ay(i)*vx(i));
    else
        ex(i) = 0;
        ey(i) = 0;
        ez(i) = 0;
    end
    if (norm(i) < 1.03*g)&&(norm(i) > 0.97*g)
        Kp = 0.5;    
        Ki = 0.005;  
    else
        Kp = 0;    
        Ki = 0;
    end
    
    if(ex(i) ~= 0 && ey(i) ~= 0 && ez(i) ~= 0) 
        %第四步误差积分   
        exInt(i+1) = exInt(i) + ex(i) * Ki * halfT * 2;
        eyInt(i+1) = eyInt(i) + ey(i) * Ki * halfT * 2;
        ezInt(i+1) = ezInt(i) + ez(i) * Ki * halfT * 2;
    else
        exInt(i+1) = 0;	
        eyInt(i+1) = 0;
        ezInt(i+1) = 0;   
    end
    %第五步互补滤波，姿态误差补偿到角速度上，修正角速度积分漂移      可以在这儿改
    gx(i+1) = gx(i) + Kp*ex(i) + exInt(i+1);
    gy(i+1) = gy(i) + Kp*ey(i) + eyInt(i+1);
    gz(i+1) = gz(i) + Kp*ez(i) + ezInt(i+1);	
    %第六步一阶龙格库塔更新四元数
    q0(i+1) = q0(i) + (-q1(i)*gx(i+1) - q2(i)*gy(i+1) - q3(i)*gz(i+1))*halfT;
    q1(i+1) = q1(i) + (q0(i)*gx(i+1) + q2(i)*gz(i+1) - q3(i)*gy(i+1))*halfT;
    q2(i+1) = q2(i) + (q0(i)*gy(i+1) - q1(i)*gz(i+1) + q3(i)*gx(i+1))*halfT;
    q3(i+1) = q3(i) + (q0(i)*gz(i+1) + q1(i)*gy(i+1) - q2(i)*gx(i+1))*halfT;
    %第七步归一化四元素
    normq(i) = sqrt(q0(i+1)*q0(i+1) + q1(i+1)*q1(i+1) + q2(i+1)*q2(i+1) + q3(i+1)*q3(i+1));
    q0(i+1) = q0(i+1) / normq(i);
    q1(i+1) = q1(i+1) / normq(i);
    q2(i+1) = q2(i+1) / normq(i);
    q3(i+1) = q3(i+1) / normq(i);

    %转化矩阵
    ux(i) = q0(i+1)*q0(i+1)+q1(i+1)*q1(i+1)-q2(i+1)*q2(i+1)-q3(i+1)*q3(i+1);%11
    uy(i) = 2.0*(q1(i+1)*q2(i+1)-q0(i+1)*q3(i+1));                          %21
    uz(i) = 2.0*(q1(i+1)*q3(i+1)+q0(i+1)*q2(i+1));                          %31
       
    wx(i) = 2.0*(q1(i+1)*q2(i+1)+q0(i+1)*q3(i+1));                          %12
    wy(i) = q0(i+1)*q0(i+1)-q3(i+1)*q3(i+1)-q1(i+1)*q1(i+1)+q2(i+1)*q2(i+1);%22
    wz(i) = 2.0*(q2(i+1)*q3(i+1)-q0(i+1)*q1(i+1));                          %32
    
    vx(i) = 2.0*(q1(i+1)*q3(i+1) - q0(i+1)*q2(i+1));                        %13
    vy(i) = 2.0*(q0(i+1)*q1(i+1) + q2(i+1)*q3(i+1));                        %23
    vz(i) = q0(i+1)*q0(i+1)-q1(i+1)*q1(i+1)-q2(i+1)*q2(i+1)+q3(i+1)*q3(i+1);%33	
    
    attitude(i,1) = asin (vy(i));%atan2(vy(i),vz(i))* 57.3;计算结果一样
    attitude(i,2) = -atan2(vx(i), vz(i));%-asin(vx(i))* 57.3;计算结果一样
    attitude(i,3) = -atan2(wx(i),ux(i));
end
attitude(:,1) = attitude(:,1)*57.3;
attitude(:,2) = attitude(:,2)*57.3;
attitude(:,3) = attitude(:,3)*57.3 + 180;

%转化成秒
% N = leng - n_start+1;
% n = 0:N;
% t = n/fs;

end
