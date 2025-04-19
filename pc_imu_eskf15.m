clear;
clc;
close all;
%% 东北天坐标系 PC IMU 组合导航，经典15维ESKF

%% 采用4.11数据
addpath('4.11');
load("INS.mat");
load("Base.mat");
load("Pol_interlated.mat");
deg2rad = pi / 180;
rad2deg = 180 / pi;
g = 9.7970;

INS = INS(1:150000,1:13);
Base = Base(1:150000,1:15);
Pol = Pol_interlated(1:150000,1:2);

for nor = 1:15
    nan = find(isnan(Base(:,11)));
    Base(nan,2:15) = Base(nan-1,2:15);
end
time_acc_1 = fix(INS(:,1)/10000);
time_Pol_1 = fix(Pol(:,1));
time_base_1 = fix(Base(:,1)/10000);

%% time alig
[Index] = time_alig(time_acc_1,time_Pol_1,time_base_1);

[start,ind] = max([time_base_1(1,1),time_acc_1(1,1),time_Pol_1(1,1)]);
if ind ==1
    start1 = (time_base_1(1,1) - time_Pol_1(1,1));
    start2 = (time_base_1(1,1) - time_acc_1(1,1));   
    
    INS = INS(start2+1:end,1:13);
    Pol = Pol_interlated(start1+1:end,1:2);

    time_acc_1 = time_acc_1(start2+1:end);
    time_Pol_1 = time_Pol_1(start1+1:end);
elseif ind ==2
    start1 = (time_acc_1(1,1) - time_base_1(1,1));
    start2 = (time_acc_1(1,1) - time_Pol_1(1,1));  

    Base = Base(start1+1:end,1:15);
    Pol = Pol_interlated(start2+1:end,1:2);
    
    time_base_1 = time_base_1(start1+1:end);
    time_Pol_1 = time_Pol_1(start2+1:end);
else
    start1 = (time_Pol_1(1,1) - time_base_1(1,1));
    start2 = (time_Pol_1(1,1) - time_acc_1(1,1));  

    Base = Base(start1+1:end,1:15);
    INS = INS(start2+1:end,1:13);

    time_base_1 = time_base_1(start1+1:end);
    time_acc_1 = time_acc_1(start2+1:end);
end
endline = min([length(time_base_1),length(time_acc_1),length(time_Pol_1)]);

Base = Base(1:endline,1:15);
INS = INS(1:endline,1:13);
Pol = Pol_interlated(1:endline,1:2);
time_base_1 = time_base_1(1:endline);
time_acc_1 = time_acc_1(1:endline);
time_Pol_1 = time_Pol_1(1:endline);

%% extract data
X_ACCL = INS(:,8) / g;%单位:g
Y_ACCL = INS(:,9) / g;
Z_ACCL = INS(:,10) / g;
acc = [X_ACCL, Y_ACCL, Z_ACCL];
X_GYRO = INS(:,5);
Y_GYRO = INS(:,6);
Z_GYRO = INS(:,7);
gyro = [X_GYRO, Y_GYRO, Z_GYRO];
X_MAG = INS(:,11)/100000;
Y_MAG = INS(:,12)/100000;
Z_MAG = INS(:,13)/100000;
pitchEKF = INS(:,3); 
rollEKF = INS(:,4);
INS_yaw = INS(:,2)';
INS_yaw(INS_yaw>=360) = INS_yaw(INS_yaw>=360)-360;
INS_yaw(INS_yaw<0) = INS_yaw(INS_yaw<0)+360;
base_heading = Base(:,11)';
pc_heading = Pol(:,2)';
pc_heading(pc_heading>=360) = pc_heading(pc_heading>=360)-360;
pc_heading(pc_heading<0) = pc_heading(pc_heading<0)+360;

u = [acc'; gyr'];
PC = [zeros(1,length(Pol));zeros(1,length(Pol));(INS_yaw-pc_heading)];
PC_time = time_Pol_1;
imu_t = time_acc_1;


%% load settings
settings = pc_imu_parameters_settings();

%% Run the GNSS-aided INS
disp('Runs the PC-aided INS');

proc_div = 0; %过程递推分频器
x = init_navigation_state(u, settings);

% 初始零偏
delta_u = zeros(6, 1);

% 初始化滤波器
[P, Q, ~, ~] = init_filter(settings);

n = length(u);

ctr_PC_data = 1;

for k=2:n
    dt = imu_t(k)-imu_t(k-1);
  
    u_h = u(:,k) - delta_u;
    
    [x(1:3), x(4:6), x(7:10)] = ch_nav_equ_local_tan(x(1:3), x(4:6), x(7:10), u_h(1:3), u_h(4:6), dt, settings.gravity);
    
    proc_div = proc_div+1;
    if proc_div == 10
        
        [F, G] = state_space_model(x, u_h, dt*proc_div);
        
        P = F*P*F' + G*Q*G';
        
        proc_div = 0;
    end

    if abs(imu_t(k) - PC_time(ctr_PC_data)) < 0.01
            
        z = PC(:, ctr_PC_data);
        M =  [0 0 0;
              0 0 0;
              0 0 1];
        H = [M zeros(3,12)];
        R = [settings.sigma_heading^2*M];
        
        % kf
        [x_kf,P_kf] = kf(x,P,F,H,z,Q,R);
        % akf
        % parameter choice
        lamda = 0.4;
        [x_akf,P_akf] = akf(x,P,F,H,z,Q,R,lamda);
        % israkf
        % parameter choice
        alfa = [0.5;zeros(length(n),1)];
        eta1 = 2; eta2 = 0.7;
        [x_israkf,P_israkf] = israkf(x,P,F,H,z,Q,R,alfa,eta1,eta2,k);
        % vbakf
        % parameter choice
        N = 10; v1 = 5; v2 = 5;
        [x_vbakf,P_vbakf] = vbakf(x,P,F,H,z,Q,R,N,v1,v2);
        % vbisrakf
        % parameter choice
        alfa = [0.5;zeros(length(n),1)];
        eta1 = 2; eta2 = 0.7; N = 10; u = 5; U = 5;
        [x_vbisrakf,P_vbisrakf,ukk,Ukk,D_R] = vbisrakf(x,P,F,H,z,Q,R,uk1k,Uk1k,gama,N,alfa,eta1,eta2,k);

        ctr_PC_data = min(ctr_PC_data+1, length(PC_time));
    end
    
    % Save the data to the output data structure
    outdata.kf = x_kf;
    outdata.akf = x_akf;
    outdata.israkf = x_israkf;
    outdata.vbakf = x_vbakf;
    outdata.vbisrakf = x_vbisrakf;
end

%% 误差计算
imuerror = base_heading - INS_yaw;
imuerror(imuerror<-180) = imuerror(imuerror<-180)+360;
imuerror(imuerror>180) = imuerror(imuerror>180)-360;

Polerror = base_heading - pc_heading;
Polerror(Polerror<-180) = Polerror(Polerror<-180)+360;
Polerror(Polerror>180) = Polerror(Polerror>180)-360;

proposederror = base_heading - outdata.vbisrakf;
proposederror(proposederror<-180) = proposederror(proposederror<-180)+360;
proposederror(proposederror>180) = proposederror(proposederror>180)-360;

AKF_error = base_heading - outdata.akf;
AKF_error(AKF_error<-180) = AKF_error(AKF_error<-180)+360;
AKF_error(AKF_error>180) = AKF_error(AKF_error>180)-360;

VBAKF_error = base_heading - outdata.vbakf;
VBAKF_error(VBAKF_error<-180) = VBAKF_error(VBAKF_error<-180)+360;
VBAKF_error(VBAKF_error>180) = VBAKF_error(VBAKF_error>180)-360;

ISRAKF_error = base_heading - outdata.israkf;
ISRAKF_error(ISRAKF_error<-180) = ISRAKF_error(ISRAKF_error<-180)+360;
ISRAKF_error(ISRAKF_error>180) = ISRAKF_error(ISRAKF_error>180)-360;

imu_max = max(imuerror);
imu_min = min(imuerror);
imu_mean = mean(imuerror);
imu_error = std(imuerror);
imu_rms = sqrt(mean(imuerror.^2));

Pol_max = max(Polerror);
Pol_min = min(Polerror);
Pol_mean = mean(Polerror);
Pol_error = std(Polerror);
Pol_rms = sqrt(mean(Polerror.^2));

proposed_max = max(proposederror);
proposed_min = min(proposederror);
proposed_mean = mean(proposederror);
proposed_std = std(proposederror);
proposed_rms = sqrt(mean(proposederror.^2));

AKF_max = max(AKF_error);
AKF_min = min(AKF_error);
AKF_mean = mean(AKF_error);
AKF_std = std(AKF_error);
AKF_rms = sqrt(mean(AKF_error.^2));

VBAKF_max = max(VBAKF_error);
VBAKF_min = min(VBAKF_error);
VBAKF_mean = mean(VBAKF_error);
VBAKF_std = std(VBAKF_error);
VBAKF_rms = sqrt(mean(VBAKF_error.^2));

ISRAKF_max = max(ISRAKF_error);
ISRAKF_min = min(ISRAKF_error);
ISRAKF_mean = mean(ISRAKF_error);
ISRAKF_std = std(ISRAKF_error);
ISRAKF_rms = sqrt(mean(ISRAKF_error.^2));

%% 绘图

figure('Name','yaw Data');
hold on;
plot(time,base_heading,'r','Linewidth',3.0);%基准
plot(time,INS_yaw,'g','Linewidth',2.0);%INS
plot(time,pc_heading,'b','Linewidth',2.0);%PC
plot(time,outdata.akf,'m','Linewidth',2.0);%AKF
plot(time,outdata.israkf,'k','Linewidth',2.0);%ISRAKF
plot(time,outdata.vbakf,'Color',[0.5,0,0.5],'Linewidth',2.0);%VBAKF
plot(time,outdata.vbisrakf,'c','Linewidth',3.0);%proposed
grid
legend('Reference','INS','PC','AKF','ISRAKF','VBAKF','Proposed');
xlabel('Time(s)');
ylabel('Heading(deg)');
hold off;

figure('Name','yaw error');
hold on;
plot(time,imuerror,'g','Linewidth',2.0);
plot(time,Polerror,'b','Linewidth',2.0);
plot(time,AKF_error,'m','Linewidth',2.0);
plot(time,ISRAKF_error,'k','Linewidth',2.0);
plot(time,VBAKF_error,'Color',[0.5,0,0.5],'Linewidth',2.0);
plot(time,proposederror,'c','Linewidth',3.0);
grid
legend('INS','PC','AKF','ISRAKF','VBAKF','Proposed');
xlabel('Time(s)');
ylabel('Heading error(deg)');
title(['航向角误差','INS',num2str(imu_rms),'PC',num2str(Pol_rms),'AKF',num2str(AKF_rms),'ISRAKF',num2str(ISRAKF_rms),'VBAKF',num2str(VBAKF_rms),'Proposed',num2str(proposed_rms)]);
hold off