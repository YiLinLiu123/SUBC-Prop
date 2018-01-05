
%% Generating Guide Curves
clear; 
clf;
points = 100; % number of data points (x,y,z) 

% y is the radially outward direcion  (along the span of the prop)
% the guidance curve of the propellers
startY = -40; % starting y  (coordinate) value for getting the correct curvture (mm) 
 % of the elipse (mm)
endY = 60; % ending Y-coordinate (for gui
y = linspace(startY,endY,points); %% radius direction, generating the list of y-coordinates for each point

% guidance curve equations
x_Upper  = ((60^2 - y.^2)./10).^(1/2); % 6^2 = x^2+10y^2
x_Lower = ((60^2 - y.^2)./30).^(1/2)* -1; % 6^2 = x^2+30y^21/4 

% translating to become the 1/4 chord length (an average chord length is
% taken) 
chordLength = x_Upper-x_Lower; %chord length
chord = x_Lower + 3/4*chordLength; %location of 1/4 chord above the y-axis 

%avg_Chord = mean(chord(1:90)); %average chord position (only took average before tampering off)
%x_Upper = x_Upper - avg_Chord; % shifted upper guidance curve
%x_Lower = x_Lower - avg_Chord; % shifted lower guidance curve 



%% Maniuplating the guidance curve for twist and skew

theta = ones(1,length(y))*60-log(y-startY+1); % test function for twist (angles in degrees)
skewAngle = 0.01 * (y-startY).^(1.5); %function for skew angle along radius (degrees)
skewDistance = tand(skewAngle) .* (y-startY); % calculating the skew distance along y-axis (mm)

%airfoil chord is in the x-z plane, with chord along the x-axis initally
% place holder vectors
x_Upper_Rotated = zeros(1,points);
x_Lower_Rotated = zeros(1,points);
z_Rotated_Upper = zeros(1,points); 
z_Rotated_Lower = zeros(1,points);


% transforming the points
for n = 1: points
    % the manuipulation is done for each x-z point
    target = [x_Upper(n)-chord(n),x_Lower(n)-chord(n);0,0]; %translate to new origin
    transformed = rotate2D(target,theta(n)); % rotation 
    x_Upper_Rotated(n) = transformed(1,1)- skewDistance(n)+chord(n); % adding skew and moving origin back
    x_Lower_Rotated(n) = transformed(1,2)- skewDistance(n)+chord(n); 
    
    z_Rotated_Upper(n) = transformed(2,1);
    z_Rotated_Lower(n) = transformed(2,2);
end

%% writing out the curve files (txt files)
% they will be saved at which directory specified here.
cd 'P:\University\DuringUni\SubC\Fluids\props\RoughProps4\AirFoils\betterProfiles' 
% location of the current folder where this code will output files to
twisted_Upper = vertcat(x_Upper_Rotated,y-startY, z_Rotated_Upper); 
twisted_Lower = vertcat(x_Lower_Rotated,y-startY,z_Rotated_Lower);
dlmwrite('twisted_Upper.txt',transpose(twisted_Upper),'delimiter','\t','precision',5);
dlmwrite('twisted_Lower.txt',transpose(twisted_Lower),'delimiter','\t','precision',5);


% plotting resultant 3d curves
axis equal
plot3(x_Lower_Rotated,y+40,z_Rotated_Lower);
hold on; 
plot3(x_Upper_Rotated,y+40,z_Rotated_Upper);

% runs the profile generation code for the beginning and end profiles
 run('ProfileExperiment.m');

%% testing

testing_One = [x_Upper(1),x_Lower(1);0,0];
 result = rotate2D(testing_One,50);
 
 
 % testing initial translated guidance curves
%{
 %test1= x_Upper - avg_Chord; 
%test2 = x_Lower - avg_Chord;

plot(x_Upper,y,'LineWidth', 2.5);
hold on
plot(x_Lower,y,'LineWidth', 2.5) ;
plot(chord,y)
plot(avg_Chord,y);
%plot(y,test1)
%plot(y,test2)
%} 
axis equal 
hold off


% verify results:


%{
test_z = [x_Upper(3), x_Lower(3)];
test_z2= [x_Upper_Rotated(3), x_Lower_Rotated(3)];
test_x = [0,0]; 
test_x2 = [x_Upper_Rotated(3),x_Lower_Rotated(3)]; 
plot(test_z,test_x)
hold on 
plot(test_z2,test_x2)
%} 


%% rotation matrix
% rotation measured from x-axis and goes ccw

function rotatedMatrix = rotate2D(input, angle)
rotation = [cosd(angle), -sind(angle); sind(angle), cosd(angle)];
rotatedMatrix = rotation * input;
end


