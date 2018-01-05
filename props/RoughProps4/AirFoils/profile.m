
filename = 'masterNormalized.txt'; %% file name for the ordinates of the airfoil 
% assume that the aifoil is normalized
delimiterIn = '\t';
normalized = importdata(filename,delimiterIn);

%% transformations 
shifted = normalized;
shifted(:,1) = -shifted(:,1); % flipping the airfoil about y-axis

% Generating the first profile (at the "origin")
first_y = 1; %point number please, location of first profile (mm along radius)
last_y = 99; % point number please

 % chord lengths for the profiles 
 first_Chord = x_Upper(first_y) - x_Lower(first_y);
 last_Chord = x_Upper(last_y) - x_Lower(last_y); 
    
 % rotation angles
 first_Theta = theta(first_y);
 last_Theta = theta(last_y);
  
% now making the transformations
    % transformation in terms of scaling the chord and then rotation
    first_Profile = shifted * first_Chord+ x_Upper(first_y); % scaling and moving to match originial profile
    first_Profile(:,1) = first_Profile(:,1)-chord(first_y); % moving the origin airfoil to 1/4 chord
    first_Profile = rotate2D(transpose(first_Profile),first_Theta); %rotation
    first_Profile = transpose(first_Profile); % geting it to be in x,y as columns
    first_Profile(:,1) = first_Profile(:,1)+chord(first_y); %moving origin back
    
    last_Profile = shifted*last_Chord;
    last_Profile(:,1) = last_Profile(:,1)+x_Upper_Rotated(last_y); 
    last_Profile = rotate2D(transpose(last_Profile),last_Theta);
    last_Profile = transpose(last_Profile);
    last_Profile(:,1) = last_Profile(:,1)+ x_Upper_Rotated(first_y);
    
    % generating the corresponding "y" components
    first_Profile_y = zeros(1, length(first_Profile))+y(first_y)-startY; 
    last_Profile_y =zeros(1, length(last_Profile))+y(last_y)-startY; 
    
    % modifying the x-compoenent for skew:
    first_Profile(:,1) = first_Profile(:,1)-skewDistance(first_y);
    last_Profile(:,1) = last_Profile(:,1) - skewDistance(last_y);
 
  % writing to a file
 first_Coordinates = horzcat(first_Profile(:,1), transpose(first_Profile_y), first_Profile(:,2));
 last_Coordinates = horzcat(last_Profile(:,1), transpose(last_Profile_y), last_Profile(:,2));
 
 
 %% writing text file
 cd 'P:\University\DuringUni\SubC\Fluids\props\RoughProps4\AirFoils\betterProfiles'
 dlmwrite('first_Profile.txt',first_Coordinates, 'delimiter','\t','precision',5);
 dlmwrite('last_Profile.txt',last_Coordinates, 'delimiter','\t','precision',5);
    
%% testing  

clf
hold on
axis equal
%{


%}

% plotting the orginial airfoil the shifted airfoil
plot(normalized(:,1), normalized(:,2));
hold on
plot(shifted(:,1),shifted(:,2));
plot(first_Profile(:,1), first_Profile(:,2)); 

%chord lines
first_Test = [twisted_Upper(1,1), twisted_Lower(1,1); twisted_Upper(3,1),twisted_Lower(3,1)];
plot(first_Test(1,:), first_Test(2,:) ); 

test_angle = atand(twisted_Upper(3,1)/twisted_Upper(1,1));
indice = find(first_Profile(2,:) == max(first_Profile(2,:))); 
test_angle_two = atand(first_Profile(2,indice)/ first_Profile(1,indice));

%{

%} 

%% rotation function (ccw)
function rotatedMatrix = rotate2D(input, angle)
rotation = [cosd(angle), -sind(angle); sind(angle), cosd(angle)];
rotatedMatrix = rotation * input;
end