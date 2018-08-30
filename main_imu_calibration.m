close all;clear all;clc;
addpath('geometry_library');
%% IMU (accelerometer & gyroscope) calibration algorithm
% Paper : D. Tedaldi et. al., "A Robust and Easy to Implement Method for IMU
% Calibration without External Equipments," IEEE ICRA, Hong Kong, 2014.
% Code author : Changhyeon Kim
% Affiliation : Intelligent Control Systems Laboratory, Seoul National University, Seoul, South Korea.
% Date        : 2018-08-30

%%
% Parameters
global freq sig_a sig_w grav T_wait T_init
freq  = 100;     % [Hz]
sig_a = 1*1e-3;  % noise density of acclerometer.
sig_w = 1*1e-3;  % noise density of gyroscope.
grav  = 9.79871; % magnitude of gravitational acceleration [m/s^2] at SNU.
% Actual nominal value : 9.80665

T_wait = 2;    % T_wait (refered in the paper), [s];
T_init = 4;      % T_init (refered in the paper), [s];

% Import file ( All data must have a column vector form. )
data = importdata('imu_good3.txt');
% data = importdata('../../save_exp/imu.txt');

data.data(1,:) = [];

t_imu = data.data(:,1).';       t_imu = t_imu - t_imu(1); acc_imu = data.data(:,2:4).';
w_imu = data.data(:,5:7).';   mag_imu = data.data(:,8:10).';
q_imu = data.data(:,11:14).'; q_imu = q_imu([2,3,4,1],:);

q_imu_init = q_imu(:,1);

% Derotate the IMU intrinsic AHRS output.
for i = 1:length(q_imu)
    q_imu(:,i) = quat_prod_kch(q_imu(:,i), quat_inv_kch(q_imu_init));
end

% Plot the time history
figure();
plot(t_imu(2:end)-t_imu(1:end-1),'.'); xlabel('seq'); ylabel('time [s]'); title('frequency history');
fprintf('time mean : %0.3f [ms], std : %0.3f [ms]\n',mean(t_imu(2:end)-t_imu(1:end-1))*1e3, std(t_imu(2:end)-t_imu(1:end-1))*1e3);

%% (1) Gyro bias correction
num_avg = floor(T_init*freq);
bias_g = [mean(w_imu(1,1:num_avg)),mean(w_imu(2,1:num_avg)),mean(w_imu(3,1:num_avg))].';

w_unbiased = w_imu-bias_g;
fprintf('Gyro bias b_g : %0.3f [rad/s], %0.3f [rad/s], %0.3f [rad/s]\n',bias_g(1), bias_g(2), bias_g(3));

figure();
subplot(3,1,1); plot(t_imu, w_unbiased(1,:)); hold on; line([1,1]*T_init,[-1,1]*100,'color',[1,0,1],'linewidth',1.5); ylabel('w_x [rad/s]'); ylim([-1,1]*5); legend('data','T_i_n_i_t'); title('Gyro bias corrected');
subplot(3,1,2); plot(t_imu, w_unbiased(2,:)); hold on; line([1,1]*T_init,[-1,1]*100,'color',[1,0,1],'linewidth',1.5);ylabel('w_y [rad/s]'); ylim([-1,1]*5);
subplot(3,1,3); plot(t_imu, w_unbiased(3,:)); hold on; line([1,1]*T_init,[-1,1]*100,'color',[1,0,1],'linewidth',1.5);ylabel('w_z [rad/s]'); xlabel('time [s]'); ylim([-1,1]*5);

%% (2) Find distinct stationary stages by using Allan variance of acceleration.
num_allan_half = floor(T_wait*freq/2);
num_allan = 2*num_allan_half + 1;
len_total = length(t_imu);

for i = 1:len_total
    if(i < num_allan_half + 1)
        var_x = var( acc_imu(1,1:i+num_allan_half) );
        var_y = var( acc_imu(2,1:i+num_allan_half) );
        var_z = var( acc_imu(3,1:i+num_allan_half) );
        
        allan_var(i) =  sqrt( var_x^2 + var_y^2 + var_z^2);
    elseif(i < len_total - num_allan_half)
        var_x = var( acc_imu(1,i-num_allan_half:i+num_allan_half) );
        var_y = var( acc_imu(2,i-num_allan_half:i+num_allan_half) );
        var_z = var( acc_imu(3,i-num_allan_half:i+num_allan_half) );
        
        allan_var(i) = sqrt( var_x^2 + var_y^2 + var_z^2);
    else
        allan_var(i) = allan_var(len_total - num_allan_half-1);
    end
    if(rem(i,floor(len_total/10)) == 1)
        fprintf('  Allan variance calculations : %0.2f [percents]\n',i/len_total*100);
    elseif(i==len_total)
        fprintf('  DONE, Allan variance calculations : %0.2f [percents]\n',i/len_total*100);
    end
end

% Difference of Allan variance to find the rising edge of the Allan
% variance time history.
for i=1:len_total
    if(i == len_total)
        diff_allan_var(i) = diff_allan_var(i-1);
    else
        diff_allan_var(i) = (allan_var(i+1) - allan_var(i));%/(t_imu(i+1) - t_imu(i));
    end
end

% Find flat (stationary) regions.
rising_flag  = false;
count_stage    = 1;      % indicator for noting an each stage.
allan_thres  = 0.1;    % threshold value for Allan variance.
transit_start  = [];     % storage for saving the start index of each stage.
transit_end    = [];     % storage for saving the end index of each stage.
stage_vec = zeros(len_total,1);

for i=1:len_total
    allan_tmp = allan_var(i);
    if(rising_flag == false && allan_tmp >= allan_thres)
        rising_flag = true;
        count_stage = count_stage + 1;
        transit_start = [transit_start, i];
    end
    if(rising_flag == true && allan_tmp < allan_thres) %transition part end
        rising_flag = false;
        transit_end   = [transit_end, i];
    end
    if(rising_flag == true)
        stage_vec(i)= count_stage - 0.5;
    else
        stage_vec(i)= count_stage;
    end
end
figure();
plot(t_imu,stage_vec,'k+');ylim([0,max(stage_vec)+1]);
hold on;
line([t_imu(transit_start);t_imu(transit_start)],repmat([0;1]*1e3,1,length(transit_start)),'Color',[1,0,1],'linewidth',1.0);
line([t_imu(transit_end);t_imu(transit_end)],repmat([0;1]*1e3,1,length(transit_end)),'Color',[0,0,1],'linewidth',1.0);

title('state history'); xlabel('time [s]'); ylabel('stage');

figure();
plot(t_imu, allan_var,'k','linewidth',2);
hold on;
line([t_imu(transit_start);t_imu(transit_start)],repmat([0;1]*1e3,1,length(transit_start)),'Color',[1,0,1],'linewidth',1.0);
line([t_imu(transit_end);t_imu(transit_end)],repmat([0;1]*1e3,1,length(transit_end)),'Color',[0,0,1],'linewidth',1.0);
ylim([0,100]);
xlabel('time [s]'); ylabel('Allan var. [m/s^2]'); title('Allan var. history of acc.');

% Extract each stage.
num_stage = floor(stage_vec(end));

acc_stage = cell(num_stage,1);
t_stage   = cell(num_stage,1);
w_trans   = cell(num_stage-1,1);
t_trans   = cell(num_stage-1,1);

for i=1:len_total
    if( rem(stage_vec(i),1) == 0) % extract stationary parts
        acc_stage{stage_vec(i)} = [acc_stage{stage_vec(i)}, acc_imu(:,i)];
        t_stage{stage_vec(i)}   = [t_stage{stage_vec(i)}, t_imu(i)];
    else % extract transition parts
        w_trans{stage_vec(i)-0.5} = [w_trans{stage_vec(i)-0.5},w_imu(:,i)];
        t_trans{stage_vec(i)-0.5}   = [t_trans{stage_vec(i)-0.5}, t_imu(i)];
    end
end

fprintf('the number of stages: %d\n',num_stage);

%% (3) Acceleration scale, misalignment and bias correction
% Average accelerometer values on each stationary stage.
% TODO: fault (movement) detector.
for i=1:num_stage
    acc_avg(:,i) = mean(acc_stage{i}.').';
end

% Iteratively calculate acc. parameter vector.
theta = zeros(9,1);
theta(:,1) = [0,0,0,1,1,1,0,0,0].';
theta_save = [];
residual_save = [];

iter_count = 1;
curr_residual = 0;
prev_residual = 1e99;
while(1)
    % Calculate residual (acc).
    r_acc = calc_residual_acc(acc_avg,theta);
    
    % Calculate Jacobian (acc).
    J_acc = calc_Jacobian_acc(acc_avg, theta);
    
    % Update the parameter vector, theta.
    d_theta = -0.25*(J_acc.'*J_acc)^-1*J_acc.'*r_acc;
    theta = theta + d_theta;
    
    % Save the parameter and residual histories.
    theta_save = [theta_save,theta];
    residual_save = [residual_save,norm(r_acc)];
    
    % Finish criterion
    curr_residual = norm(r_acc.^2);
    if(prev_residual - curr_residual <= 1e-10)
        break;
    end
    prev_residual = curr_residual;
    
    % Increase the counter.
    iter_count = iter_count + 1;
    %     fprintf(' ACC iteration : %d\n', iter_count);
end

a_yz = theta(1);
a_zy = theta(2);
a_zx = theta(3);
sx = theta(4);
sy = theta(5);
sz = theta(6);
bx = theta(7);
by = theta(8);
bz = theta(9);

bias_a = [bx,by,bz].';

figure();
plot(residual_save); title('residual');

fprintf('Optimized values - a_yz: %0.3f, a_zy: %0.3f, a_zx: %0.3f, sx: %0.3f, sy: %0.3f, sz: %0.3f, bx: %0.3f, by: %0.3f, bz: %0.3f\n',a_yz,a_zy,a_zx,sx,sy,sz,bx,by,bz);

Sa = diag([sx,sy,sz]);
Ta = [1,-a_yz,a_zy;...
    0,1,-a_zx;...
    0,0,1];

acc_avg_o = Ta*Sa*(acc_avg + bias_a);

figure();
plot3(acc_avg(1,:),acc_avg(2,:),acc_avg(3,:),'bo');
grid on; hold on;
plot3(acc_avg_o(1,:),acc_avg_o(2,:),acc_avg_o(3,:),'r+');
axis square;
xlim([-1,1]*11); ylim([-1,1]*11); zlim([-1,1]*11); title('3d acc. plot');

acc_o = Ta*Sa*(acc_imu + [bx,by,bz].');

figure();
subplot(3,1,1);
plot(t_imu, acc_imu(1,:),'k'); hold on;
plot(t_imu,acc_o(1,:),'r');
ylabel('acc x [m/s^2]'); ylim([-1,1]*9.81*2);
line([t_imu(transit_start);t_imu(transit_start)],repmat([-1;1]*1e3,1,length(transit_start)),'Color',[1,0,1],'linewidth',1.0);
line([t_imu(transit_end);t_imu(transit_end)],repmat([-1;1]*1e3,1,length(transit_end)),'Color',[0,0,1],'linewidth',1.0);
subplot(3,1,2);
plot(t_imu, acc_imu(2,:),'k'); hold on;
plot(t_imu,acc_o(2,:),'r');
ylabel('acc y [m/s^2]'); ylim([-1,1]*9.81*2);
line([t_imu(transit_start);t_imu(transit_start)],repmat([-1;1]*1e3,1,length(transit_start)),'Color',[1,0,1],'linewidth',1.0);
line([t_imu(transit_end);t_imu(transit_end)],repmat([-1;1]*1e3,1,length(transit_end)),'Color',[0,0,1],'linewidth',1.0);
subplot(3,1,3);
plot(t_imu, acc_imu(3,:),'k'); hold on;
plot(t_imu,acc_o(3,:),'r');
ylabel('acc z [m/s^2]'); ylim([-1,1]*9.81*2);
line([t_imu(transit_start);t_imu(transit_start)],repmat([-1;1]*1e3,1,length(transit_start)),'Color',[1,0,1],'linewidth',1.0);
line([t_imu(transit_end);t_imu(transit_end)],repmat([-1;1]*1e3,1,length(transit_end)),'Color',[0,0,1],'linewidth',1.0);


%% (4) Gyroscope scale, misalignment correction
% Nonimal vectors from calibrated accelerometer measurements.
acc_avg_o = Ta*Sa*(acc_avg + bias_a);
for i=1:num_stage
    u_ref(:,i) = acc_avg_o(:,i) / norm(acc_avg_o(:,i));
end

% Calculate rotation matrices between i-th and (i+1)-th stages by
% Select relative rotations between stages.
R_selected = cell(num_stage-1,1);
for i = 1:num_stage-1
    q_tmp = [1,0,0,0].';
    for j = transit_start(i):transit_end(i)
        dt = t_imu(j+1) - t_imu(j);
        %         w_temp = diag([-1,-1,1])*w_unbiased(:,j);
        w_tmp = w_unbiased(:,j);
        
        k1 = quat_derivative_kch(q_tmp,           w_tmp);
        k2 = quat_derivative_kch(q_tmp + k1*dt/2, w_tmp);
        k3 = quat_derivative_kch(q_tmp + k2*dt/2, w_tmp);
        k4 = quat_derivative_kch(q_tmp + k3*dt,   w_tmp);
        
        q_tmp = q_tmp + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    end
    
    % Normalize the quaternion.
    q_tmp = q_tmp/norm(q_tmp);
    R_selected{i} = quat_to_dcm_kch(quat_inv_kch(q_tmp));
end

% Rotate nominal vector
u_integral = zeros(3,num_stage-1);
u_diff     = zeros(1,num_stage-1);
for i=1:num_stage-1
    if(i==1)
        u_integral(:,i) = u_ref(:,i);
        u_integral(:,i+1) = R_selected{i}*u_ref(:,i);
    else
        u_integral(:,i+1) = R_selected{i}*u_ref(:,i);
    end
    u_diff(i) = norm(u_integral(:,i)-u_ref(:,i));
end

figure();
plot(u_diff);


%% (5) RK4 integration gyro ???

% q_RK = zeros(4,length(t_imu));
%
% norm_q = zeros(1,length(t_imu));
% norm_q(1)=1;
% q_RK(:,1) = [1;0;0;0];
%
% % Runge-Kutta 4th order integration
% for i=1:length(t_imu)-1
%    dt = t_imu(i+1) - t_imu(i);
%    w_tmp = diag([-1,-1,1])*w_unbiased(:,i); % WHY ?????????????????
%
%
%    k1 = quat_derivative_kch(q_RK(:,i), w_tmp);
%    k2 = quat_derivative_kch(q_RK(:,i)+k1*dt/2, w_tmp);
%    k3 = quat_derivative_kch(q_RK(:,i)+k2*dt/2, w_tmp);
%    k4 = quat_derivative_kch(q_RK(:,i)+k3*dt,   w_tmp);
%
%    q_RK(:,i+1) = q_RK(:,i) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
%
%    q_RK(:,i+1) = q_RK(:,i+1)/norm(q_RK(:,i+1));
% end
%
%
% figure();
% subplot(4,1,1); plot(t_imu, q_imu(1,:),'k'); hold on; plot(t_imu, q_RK(1,:),'r'); grid on; grid minor; legend('truth','RK4'); ylim([-1.5,1.5]); title('quaternion dead-reckoning (Runge-Kutta 4th)');
% subplot(4,1,2); plot(t_imu, q_imu(2,:),'k'); hold on; plot(t_imu, q_RK(2,:),'r'); grid on; grid minor; ylim([-1.5,1.5]);
% subplot(4,1,3); plot(t_imu, q_imu(3,:),'k'); hold on; plot(t_imu, q_RK(3,:),'r'); grid on; grid minor; ylim([-1.5,1.5]);
% subplot(4,1,4); plot(t_imu, q_imu(4,:),'k'); hold on; plot(t_imu, q_RK(4,:),'r'); grid on; grid minor; ylim([-1.5,1.5]);