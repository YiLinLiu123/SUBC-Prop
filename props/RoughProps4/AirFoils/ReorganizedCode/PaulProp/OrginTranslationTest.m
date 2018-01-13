%% simple shifting origin test

% right triangle of side lengh 1 with corner at (1,1)
% x,y columns
clf
clear
triangle = [1,1;5,1;1,3;1,1];
plot(triangle(:,1), triangle(:,2));
hold on; 

% moving origin
triangle(:,1)= triangle(:,1) - 1;
triangle(:,2) = triangle(:,2)- 1;

% rotate 90 degrees:
triangle_Rotated = transpose(rotate2D(transpose(triangle),90));

%shift origin back:
triangle_Rotated(:,1)= triangle_Rotated(:,1)+1;
triangle_Rotated(:,2)= triangle_Rotated(:,2)+1;

plot(triangle_Rotated(:,1),triangle_Rotated(:,2));
xlim([-10,10]);
ylim([-10,10]);


axis equal 
%% rotation matrix
% rotation measured from x-axis and goes ccw

function rotatedMatrix = rotate2D(input, angle)
rotation = [cosd(angle), -sind(angle); sind(angle), cosd(angle)];
rotatedMatrix = rotation *input;
end