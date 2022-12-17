%% ME C231A Project - Han Nguyen Fall 2022
%% Capture Data for Movement
data = csvread('move_from_zero.csv');
joint_angles = data(:,1:7);
joint_velocities = data(:,8:14);
endpoint_position = data(:,15:21);
endpoint_velocity = data(:,22:27);
l0 = data(:,28:34);
l1 = data(:,35:41);
l2 = data(:,42:48);
l3 = data(:,49:55);
l4 = data(:,56:62);
l5 = data(:,63:69);
l6 = data(:,70:76);
rh = data(:,77:83);
time = data(:,84);

%% Joint Angles over Time
for i=1:7
    subplot(4,2,i);
    plot(time,joint_angles(:,i))
    title(sprintf('Joint %d vs t', i))
end

%% Test Forecast Using Built-in Endpoint Velocity
endpoint_forecast = zeros(5098,3);
endpoint_forecast(1,:) = endpoint_position(1,1:3);
for i=1:length(endpoint_position)-1
    dT = time(i+1) - time(i);
    endpoint_forecast(i+1,:) = endpoint_forecast(i,:) + endpoint_velocity(i,1:3)*dT;
end

scatter3(endpoint_position(:,1),endpoint_position(:,2),endpoint_position(:,3));
hold on
plot3(endpoint_position(1,1),endpoint_position(1,2),endpoint_position(1,3),'x','LineWidth',10);
plot3(endpoint_position(end,1),endpoint_position(end,2),endpoint_position(end,3),'x','LineWidth',10);
scatter3(endpoint_forecast(:,1),endpoint_forecast(:,2),endpoint_forecast(:,3));
plot3(endpoint_forecast(1,1),endpoint_forecast(1,2),endpoint_forecast(1,3),'diamond','LineWidth',10);
plot3(endpoint_forecast(end,1),endpoint_forecast(end,2),endpoint_forecast(end,3),'diamond','LineWidth',10);
hold off
title('Trajectory 3D Comparison of Data vs. Forecast with Given End-Effector Velocity');
legend('Data','Data Start', 'Data End','Forecast', 'Forecast Start', 'Forecast End');

%% Test with Spatial Jacobian for Each Matrix
twists = zeros(6,7);
joints = [l0 l1 l2 l3 l4 l5 l6];
J = zeros(6,7);

endpoint_forecast = zeros(5098,3);
endpoint_forecast(1,:) = endpoint_position(1,1:3);

for i=1:length(endpoint_position)-1
    dT = time(i+1) - time(i);
    theta = joint_angles(i,:);
    theta_dot = joint_velocities(i,:);
    for j=1:7
        joint = joints(i,1+(j-1)*7:1+(j-1)*7+6);
        q = joint(1:3);
        rotation_matrix = quat2rotm(joint(4:7));
        w = rotation_matrix(:,3);
        v = -cross(w,q);
        twists(1:3,j) = v;
        twists(4:6,j) = w;
        Ad = eye(6);
        for x=1:j-1
            Ad = Ad * adjoint(twists(:,x),theta(x));
        end
        xi_prime = Ad*twists(:,j);
        J(:,j) = xi_prime;
    end
    Vs = J*theta_dot';
    endpoint_forecast(i+1,:) = endpoint_forecast(i,:) + Vs(1:3)'*dT;
end

scatter3(endpoint_position(:,1),endpoint_position(:,2),endpoint_position(:,3));
hold on
plot3(endpoint_position(1,1),endpoint_position(1,2),endpoint_position(1,3),'x','LineWidth',10);
plot3(endpoint_position(end,1),endpoint_position(end,2),endpoint_position(end,3),'x','LineWidth',10);
scatter3(endpoint_forecast(:,1),endpoint_forecast(:,2),endpoint_forecast(:,3));
plot3(endpoint_forecast(1,1),endpoint_forecast(1,2),endpoint_forecast(1,3),'diamond','LineWidth',10);
plot3(endpoint_forecast(end,1),endpoint_forecast(end,2),endpoint_forecast(end,3),'diamond','LineWidth',10);
hold off
title('Trajectory 3D Comparison of Data vs. Forecast with Spatial  Jacobian');
legend('Data','Data Start', 'Data End','Forecast', 'Forecast Start', 'Forecast End');

%% Test Forward Kinematics Forecast
gst0 = [0.0016  0.0054  1.      1.0161;
       -0.9851 -0.1719  0.0025  0.1594;
        0.1719 -0.9851  0.005   0.32;
        0.      0.      0.      1];
twists =  [0.0000   -0.3164   -0.0003   -0.3168   -0.0010   -0.3158    0.0022;
    0.0003    0.0001    0.3167    0.0001    0.3206    0.0051    0.3079;
   -0.0000    0.0805   -0.1928    0.4817   -0.0235    0.8813   -0.1654;
    0.0008    0.0001    1.0000   -0.0017    1.0000   -0.0018    1.0000;
    0.0014    1.0000    0.0010    1.0000    0.0028    1.0000   -0.0037;
    1.0000   -0.0002    0.0001   -0.0014   -0.0055   -0.0064    0.0067];

endpoint_forecast = zeros(5098,3);
endpoint_forecast(1,:) = endpoint_position(1,1:3);

for i=1:length(endpoint_position)
    theta = joint_angles(i,:);
    gst = prod_exp(twists, theta) * gst0;
    pos = gst(1:3,4)';
    endpoint_forecast(i,:) = pos';
end

scatter3(endpoint_position(:,1),endpoint_position(:,2),endpoint_position(:,3));
hold on
plot3(endpoint_position(1,1),endpoint_position(1,2),endpoint_position(1,3),'x','LineWidth',10);
plot3(endpoint_position(end,1),endpoint_position(end,2),endpoint_position(end,3),'x','LineWidth',10);
scatter3(endpoint_forecast(:,1),endpoint_forecast(:,2),endpoint_forecast(:,3));
plot3(endpoint_forecast(1,1),endpoint_forecast(1,2),endpoint_forecast(1,3),'diamond','LineWidth',10);
plot3(endpoint_forecast(end,1),endpoint_forecast(end,2),endpoint_forecast(end,3),'diamond','LineWidth',10);
hold off
title('Trajectory 3D Comparison of Data vs. Forecast with Spatial  Jacobian');
legend('Data','Data Start', 'Data End','Forecast', 'Forecast Start', 'Forecast End');
%% MPC xOpt Data Trial 1: Ts = 0.1, N = 10, Q = np.eye(3)*100
close all; clear; clc;
Trial = 1;
Ts = 0.1;
N = 10;
Q = 'np.eye(3)*100';
gst0 = [0.0016  0.0054  1.      1.0161;
       -0.9851 -0.1719  0.0025  0.1594;
        0.1719 -0.9851  0.005   0.32;
        0.      0.      0.      1];
twists =  [0.0000   -0.3164   -0.0003   -0.3168   -0.0010   -0.3158    0.0022;
    0.0003    0.0001    0.3167    0.0001    0.3206    0.0051    0.3079;
   -0.0000    0.0805   -0.1928    0.4817   -0.0235    0.8813   -0.1654;
    0.0008    0.0001    1.0000   -0.0017    1.0000   -0.0018    1.0000;
    0.0014    1.0000    0.0010    1.0000    0.0028    1.0000   -0.0037;
    1.0000   -0.0002    0.0001   -0.0014   -0.0055   -0.0064    0.0067];

data = csvread('mpc_data.csv');

endpoint_position = zeros(11,3);
for i=1:size(data,2)
    joint_angles = data(:, i);
    gst = prod_exp(twists, joint_angles) * gst0;
    pos = gst(1:3,4)';
    endpoint_position(i,:) = pos;
end

scatter3(endpoint_position(:,1),endpoint_position(:,2),endpoint_position(:,3));
hold on
plot3(endpoint_position(1,1),endpoint_position(1,2),endpoint_position(1,3),'x','LineWidth',10);
plot3(endpoint_position(end,1),endpoint_position(end,2),endpoint_position(end,3),'x','LineWidth',10);
plot3(0.5,0.5,0.5,'diamond','LineWidth',10);
title(sprintf('MPC Controller Trial %d: Ts = %0.1f, N = %d, Q = %s',Trial,Ts,N,Q));
legend('xOpt Computed','xOpt Start', 'xOpt End','Goal');

%% MPC xOpt Data Trial 2
close all; clear; clc;
Trial = 2;
Ts = 1;
N = 10;
Q = 'np.eye(3)*10';
U_lim = 0.1;
gst0 = [0.0016  0.0054  1.      1.0161;
       -0.9851 -0.1719  0.0025  0.1594;
        0.1719 -0.9851  0.005   0.32;
        0.      0.      0.      1];
twists =  [0.0000   -0.3164   -0.0003   -0.3168   -0.0010   -0.3158    0.0022;
    0.0003    0.0001    0.3167    0.0001    0.3206    0.0051    0.3079;
   -0.0000    0.0805   -0.1928    0.4817   -0.0235    0.8813   -0.1654;
    0.0008    0.0001    1.0000   -0.0017    1.0000   -0.0018    1.0000;
    0.0014    1.0000    0.0010    1.0000    0.0028    1.0000   -0.0037;
    1.0000   -0.0002    0.0001   -0.0014   -0.0055   -0.0064    0.0067];

x = csvread('state_input_constraints_x.csv');
u = csvread('state_input_constraints_input.csv');
y = csvread('state_input_constraints_output.csv')';

figure;
scatter3(y(:,1),y(:,2),y(:,3));
hold on
plot3(y(1,1),y(1,2),y(1,3),'x','LineWidth',10);
plot3(y(end,1),y(end,2),y(end,3),'x','LineWidth',10);
plot3(0.5,0.5,0.5,'diamond','LineWidth',10);
title(sprintf('MPC Controller Trial %d: Ts = %0.1f, N = %d, Q = %s',Trial,Ts,N,Q));
legend('yOpt Computed','yOpt Start', 'yOpt End','Goal');

figure;
time = 0:Ts:N;
% Joint Angles over Time
for i=1:7
    subplot(4,2,i);
    plot(time,x(i,:))
    title(sprintf('Joint %d vs Time', i))
end

figure;
% Joint Input over Time
for i=1:7
    subplot(4,2,i);
    plot(time,u(i,:))
    hold on
    ylim([-U_lim*2, U_lim*2]);
    plot([time(1) time(end)], [U_lim U_lim], ":", 'color', 'r');
    plot([time(1) time(end)], [-U_lim -U_lim], ":", 'color', 'r');
    title(sprintf('Joint Input %d vs Time', i))
end

%% MPC xOpt Data Trial State - Input
close all; clear; clc;
Ts = 1;
N = 10;
Q = 'np.eye(3)*10';
U_lim = 0.1;
gst0 = [0.0016  0.0054  1.      1.0161;
       -0.9851 -0.1719  0.0025  0.1594;
        0.1719 -0.9851  0.005   0.32;
        0.      0.      0.      1];
twists =  [0.0000   -0.3164   -0.0003   -0.3168   -0.0010   -0.3158    0.0022;
    0.0003    0.0001    0.3167    0.0001    0.3206    0.0051    0.3079;
   -0.0000    0.0805   -0.1928    0.4817   -0.0235    0.8813   -0.1654;
    0.0008    0.0001    1.0000   -0.0017    1.0000   -0.0018    1.0000;
    0.0014    1.0000    0.0010    1.0000    0.0028    1.0000   -0.0037;
    1.0000   -0.0002    0.0001   -0.0014   -0.0055   -0.0064    0.0067];

x = csvread('state_input_constraints_x.csv');
u = csvread('state_input_constraints_input.csv');
y = csvread('state_input_constraints_output.csv')';

figure;
scatter3(y(:,1),y(:,2),y(:,3));
hold on
plot3(y(1,1),y(1,2),y(1,3),'x','LineWidth',10);
plot3(y(end,1),y(end,2),y(end,3),'x','LineWidth',10);
plot3(0.5,0.5,0.5,'diamond','LineWidth',10);
title(sprintf('MPC Controller Trial with State and Input Constraints: Ts = %0.1f, N = %d, Q = %s',Ts,N,Q));
legend('yOpt Computed','yOpt Start', 'yOpt End','Goal');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');

figure;
time = 0:Ts:N;
% Joint Angles over Time
for i=1:7
    subplot(4,2,i);
    plot(time,x(i,:))
    title(sprintf('Joint %d vs Time', i))
end

figure;
% Joint Input over Time
for i=1:7
    subplot(4,2,i);
    plot(time,u(i,:))
    hold on
    ylim([-U_lim*2, U_lim*2]);
    plot([time(1) time(end)], [U_lim U_lim], ":", 'color', 'r');
    plot([time(1) time(end)], [-U_lim -U_lim], ":", 'color', 'r');
    title(sprintf('Joint Input %d vs Time', i))
end

%% MPC xOpt Data Trial Orientation Constraints
close all; clear; clc;
Ts = 1;
N = 10;
Q = 'np.eye(3)*100';
U_lim = 0.2;
gst0 = [0.0016  0.0054  1.      1.0161;
       -0.9851 -0.1719  0.0025  0.1594;
        0.1719 -0.9851  0.005   0.32;
        0.      0.      0.      1];
twists =  [0.0000   -0.3164   -0.0003   -0.3168   -0.0010   -0.3158    0.0022;
    0.0003    0.0001    0.3167    0.0001    0.3206    0.0051    0.3079;
   -0.0000    0.0805   -0.1928    0.4817   -0.0235    0.8813   -0.1654;
    0.0008    0.0001    1.0000   -0.0017    1.0000   -0.0018    1.0000;
    0.0014    1.0000    0.0010    1.0000    0.0028    1.0000   -0.0037;
    1.0000   -0.0002    0.0001   -0.0014   -0.0055   -0.0064    0.0067];

x = csvread('orientation_constraints_x.csv');
u = csvread('orientation_constraints_input.csv');
y = csvread('orientation_constraints_output.csv')';

figure;
scatter3(y(:,1),y(:,2),y(:,3));
hold on
plot3(y(1,1),y(1,2),y(1,3),'x','LineWidth',10);
plot3(y(end,1),y(end,2),y(end,3),'x','LineWidth',10);
plot3(0.747, 0.312, 0.77,'diamond','LineWidth',10);
title(sprintf('MPC Controller Trial with Orientation Constraints: Ts = %0.1f, N = %d, Q = %s',Ts,N,Q));
legend('yOpt Computed','yOpt Start', 'yOpt End','Goal');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');

figure;
time = 0:Ts:N;
% Joint Angles over Time
for i=1:7
    subplot(4,2,i);
    plot(time,x(i,:))
    title(sprintf('Joint %d vs Time', i))
end

figure;
% Joint Input over Time
for i=1:7
    subplot(4,2,i);
    plot(time,u(i,:))
    hold on
    ylim([-U_lim*2, U_lim*2]);
    plot([time(1) time(end)], [U_lim U_lim], ":", 'color', 'r');
    plot([time(1) time(end)], [-U_lim -U_lim], ":", 'color', 'r');
    title(sprintf('Joint Input %d vs Time', i))
end

%% Rotation Axis
x_axis = [];
y_axis = [];
z_axis = [];
for i=1:size(x,2)
    joint_angles = x(:,i);
    gst = prod_exp(twists, joint_angles);
    x_axis = [x_axis gst(1:3,1)];
    y_axis = [y_axis gst(1:3,2)];
    z_axis = [z_axis gst(1:3,3)];
end
%% MPC xOpt with Path Constraints
close all; clear; clc;
Ts = 1;
N = 10;
Q = 'np.eye(3)*100';
U_lim = 0.2;
gst0 = [0.0016  0.0054  1.      1.0161;
       -0.9851 -0.1719  0.0025  0.1594;
        0.1719 -0.9851  0.005   0.32;
        0.      0.      0.      1];
twists =  [0.0000   -0.3164   -0.0003   -0.3168   -0.0010   -0.3158    0.0022;
    0.0003    0.0001    0.3167    0.0001    0.3206    0.0051    0.3079;
   -0.0000    0.0805   -0.1928    0.4817   -0.0235    0.8813   -0.1654;
    0.0008    0.0001    1.0000   -0.0017    1.0000   -0.0018    1.0000;
    0.0014    1.0000    0.0010    1.0000    0.0028    1.0000   -0.0037;
    1.0000   -0.0002    0.0001   -0.0014   -0.0055   -0.0064    0.0067];

x = csvread('path_constraints_x.csv');
u = csvread('path_constraints_input.csv');
y = csvread('path_constraints_output.csv')';

figure;
scatter3(y(:,1),y(:,2),y(:,3));
hold on
plot3(y(1,1),y(1,2),y(1,3),'x','LineWidth',10);
plot3(y(end,1),y(end,2),y(end,3),'x','LineWidth',10);
plot3(0.6935683264807795-0.5, 0.16701999999999667, 0.5184984909676793,'diamond','LineWidth',10);
title(sprintf('MPC Controller Trial with Path Constraints: Ts = %0.1f, N = %d, Q = %s',Ts,N,Q));
legend('yOpt Computed','yOpt Start', 'yOpt End','Goal');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
zlim([0.3 0.7]);
ylim([0.1 0.2]);

figure;
time = 0:Ts:N;
% Joint Angles over Time
for i=1:7
    subplot(4,2,i);
    plot(time,x(i,:))
    title(sprintf('Joint %d vs Time', i))
end

figure;
% Joint Input over Time
for i=1:7
    subplot(4,2,i);
    plot(time,u(i,:))
    hold on
    ylim([-U_lim*2, U_lim*2]);
    plot([time(1) time(end)], [U_lim U_lim], ":", 'color', 'r');
    plot([time(1) time(end)], [-U_lim -U_lim], ":", 'color', 'r');
    title(sprintf('Joint Input %d vs Time', i))
end
%% Helper Functions
function S = skew_3d(omega)
    S = [[0, -omega(3), omega(2)]
         [omega(3), 0, -omega(1)]
         [-omega(2), omega(1),0]];
end

function R = rotation_3d(omega, theta)
    hat_u = skew_3d(omega);
    theta = theta * norm(omega);
    hat_u = hat_u / norm(omega);
    R = eye(3) + hat_u * sin(theta) + hat_u * hat_u * (1 - cos(theta));
end

function [R,p] = homog_3d(xi, theta)
    if all(size(xi) == [1 6])
        xi = xi';
    end
    v = xi(1:3);
    w = xi(4:6);
    I = eye(3);
    R = rotation_3d(w, theta);
    p = (1/norm(w)^2) * ((I-R)*skew_3d(w) * v + theta*(w*w') * v);
end

function A = adjoint(xi, theta)
    A = zeros(6,6);
    [R,p] = homog_3d(xi, theta);
    A(1:3,1:3) = R;
    A(1:3,4:6) = skew_3d(p)*R;
    A(4:6,4:6) = R;
end

function G = prod_exp(xi, theta)
    G = eye(4);
    for i=1:size(xi,2)
        xi_i = xi(:, i);
        theta_i = theta(i);
        [R, p] = homog_3d(xi_i, theta_i);
        g_i = eye(4);
        g_i(1:3, 1:3) = R;
        g_i(1:3, 4) = p;
        G = G * g_i;
    end
end

