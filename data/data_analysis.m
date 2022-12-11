%% ME C231A Project - Han Nguyen Fall 2022
%% Capture Data for Movement - Transforms here are inverted so they can't be used.
data = csvread('move_from_zero.csv');
joint_angles = data(:,1:7);
joint_velocities = data(:,8:14);
endpoint_position = data(:,15:21);
endpoint_velocity = data(:,22:27);
time = data(:,84);
%% Test Forecast Using Built-in Endpoint Velocity
dT = mean(diff(time));
endpoint_forecast = endpoint_position(1,1:3);
for i=1:length(endpoint_position)-1
    forecast = endpoint_forecast(i,:) + endpoint_velocity(i,1:3)*dT;
    endpoint_forecast = [endpoint_forecast; forecast];
end
figure;
scatter3(endpoint_position(:,1),endpoint_position(:,1),endpoint_position(:,1));
hold on
scatter3(endpoint_forecast(:,1),endpoint_forecast(:,1),endpoint_forecast(:,1));
hold off
title('Trajectory 3D Comparison of Data vs. Forecast');
legend('Data','Forecast');

figure;
plot(endpoint_position(1,:),endpoint_position(2,:));
hold on
plot(endpoint_forecast(1,:),endpoint_forecast(2,:));
hold off
title('Trajectory 2D Comparison of Data vs. Forecast');
legend('Data','Forecast');

%% Joint Angles over Time
for i=1:7
    subplot(4,2,i);
    plot(time,joint_angles(:,i))
end

%%
