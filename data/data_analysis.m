%% ME C231A Project - Han Nguyen Fall 2022
%% Capture Data for Movement - Transforms here are inverted so they can't be used.
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
%%
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
    v = xi(1:3);
    if all(size(v) == [1 3])
        v = v';
    end
    w = xi(4:6);
    I = eye(3);
    R = rotation_3d(w, theta);
    p = (1/norm(w)^2) * ((I-R)*skew_3d(w) * v + theta*(w'*w) * v);
end

function A = adjoint(xi, theta)
    A = zeros(6,6);
    [R,p] = homog_3d(xi, theta);
    A(1:3,1:3) = R;
    A(1:3,4:6) = skew_3d(p)*R;
    A(4:6,4:6) = R;
end
