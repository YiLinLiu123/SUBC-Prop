%% Set Up
clear
% files will be saved at which directory specified here.
cd 'P:\University\DuringUni\SubC\Fluids\props\RoughProps4\AirFoils\ReorganizedCode\PaulProp'

% importing airFoil
filename = 'masterNormalized.txt'; %% file name for the ordinates of the airfoil 
% assume that the aifoil is normalized
delimiterIn = '\t';
normalizedFoil = importdata(filename,delimiterIn);


% guidance Curve Equations
    % y vector (radial vector for the blade)
    startY=1; % starting y  (coordinate) value for getting the correct curvture (mm) 
     % of the elipse (mm)
    endY = 100; % ending Y-coordinate
    points = 99;
    y = linspace(startY,endY,points); %% radius direction, generating the list of y-coordinates 

offset = 70; % offset for the elipse (so that we get growing, then shrinking) 
x_Upper = 100-y;
x_Lower = y;
%{
x_Upper  = ((130^2 - (2.*y-offset).^2)./36).^(1/2);
x_Lower = ((130^2 - (2.*y-offset).^2)./100).^(1/2)* -1; 
%}


% 1/4 chord locations:
quarterChord = x_Lower + 3/4*(x_Upper-x_Lower);

% setting up the airfoil:
fliped_Foil= normalizedFoil;
fliped_Foil(:,1) = -fliped_Foil(:,1); % flipping the airfoil about y-axis


%% Testing for Airfoil
clf
 
% plotting the airfoil
title('NACA-2411 Airfoil')
plot(normalizedFoil(:,1),normalizedFoil(:,2), 'r');
hold on
axis equal
%{
plot(fliped_Foil(:,1), fliped_Foil(:,2),'b'); 
scatter(fliped_Foil(:,1), fliped_Foil(:,2),'b');
plot([-2,2],[0,0]); 
hold off
%} 

%% Test initial guide Curves
% plotting the initial guidance curves
title('Flat Propellor Curves')
plot(y,x_Upper,'r');
hold on 
plot(y,x_Lower,'b');
xlabel('y-coordinates (mm)') % x-axis label
ylabel('x-coordinates (mm)')
axis equal 


%% Transformations

% functions
twistAngle = ones(length(y)).*60; % test function for twist (angles in degrees)
%twistAngle = 10*y; 
skewAngle = 0.01 * (y-startY).^(1.5); %function for skew angle along radius (degrees)
skewDistance = zeros(length(y)); % calculating the skew distance along y-axis (mm)

% place holders for curve transformations:
x_Upper_Rotated = zeros(1,points);
x_Lower_Rotated = zeros(1,points);
z_Rotated_Upper = zeros(1,points); 
z_Rotated_Lower = zeros(1,points);

% airfoil profiles
% Generating the first profile (at the "origin")
first_y = 1; %point number please, location of first profile (mm along radius)
last_y = 98; % point number please


% master for loop tranformations
for n= 1: points
      % the manuipulation is done for each x-z point
    target = [x_Upper(n),x_Lower(n);0,0]; 
    
    % moving origin to 1/4 chord length
    target(1,1) = target(1,1) - quarterChord(n); 
    target(1,2) = target(1,2) - quarterChord(n);
    
    % rotation
    rotation = [cosd(twistAngle(n)), -sind(twistAngle(n)); sind(twistAngle(n)), cosd(twistAngle(n))];
    transformed = rotation * target; 
    
    % re-translation
    x_Upper_Rotated(n) = transformed(1,1)- skewDistance(n)+ quarterChord(n); % adding skew and re-translation
    x_Lower_Rotated(n) = transformed(1,2)- skewDistance(n)+ quarterChord(n);
    z_Rotated_Upper(n) = transformed(2,1);
    z_Rotated_Lower(n) = transformed(2,2);
    
    % airfoil profile transformations
    if((n == first_y) || (n == last_y))
        %scale oordinates
        chord_Length = x_Upper(n)-x_Lower(n);
        transformed_Foil= fliped_Foil .* chord_Length;
        
        % translating to match initial guidance curve
        transformed_Foil(:,1) = transformed_Foil(:,1) + x_Upper(n);
        
        %now moving to quarter chord
        transformed_Foil(:,1) = transformed_Foil(:,1) - quarterChord(n); 
        
        % now rotated
        transformed_Foil = transpose( rotation * transformed_Foil.'); 
        
        % now move back the "origin"
        transformed_Foil(:,1) = transformed_Foil(:,1)- skewDistance(n)+ quarterChord(n);
        
        % storing information
        if(n == first_y)
           first_Profile = transformed_Foil;  
        else 
           last_Profile = transformed_Foil;
        end
        
        
    end 
    
    
    
end 

%% writing out the curve files (txt files)

% location of the current folder where this code will output files to
twisted_Upper = vertcat(x_Upper_Rotated,y-startY, z_Rotated_Upper); 
twisted_Lower = vertcat(x_Lower_Rotated,y-startY,z_Rotated_Lower);

% generating the corresponding "y" components
first_Profile_y = zeros(1, length(first_Profile))+y(first_y)-startY; 
last_Profile_y =zeros(1, length(last_Profile))+y(last_y)-startY; 

%formatting profiles
 first_Coordinates = horzcat(first_Profile(:,1), transpose(first_Profile_y), first_Profile(:,2));
 last_Coordinates = horzcat(last_Profile(:,1), transpose(last_Profile_y), last_Profile(:,2));

 %{
 dlmwrite('twisted_Upper.txt',transpose(twisted_Upper),'delimiter','\t','precision',5);
 dlmwrite('twisted_Lower.txt',transpose(twisted_Lower),'delimiter','\t','precision',5);
 dlmwrite('first_Profile.txt',first_Coordinates, 'delimiter','\t','precision',5);
 dlmwrite('last_Profile.txt',last_Coordinates, 'delimiter','\t','precision',5);
%}

%% Testing Code for Transformes 

clf
% plotting resultant 3d curves
axis equal
plot3(x_Lower_Rotated,y-startY,z_Rotated_Lower,'b');
hold on; 
plot3(x_Upper_Rotated,y-startY,z_Rotated_Upper,'r');
plot3(first_Coordinates(:,1), first_Coordinates(:,2), first_Coordinates(:,3));
plot3(last_Coordinates(:,1), last_Coordinates(:,2), last_Coordinates(:,3));
xlabel('x-coordinates (mm)')
ylabel('y-coordinates (mm)')
zlabel('z-coordinates (mm)') 
title('all curves for propeller')
%legend('trailing edge curve','leading edge Curve','first Cross Section', 'Second Cross Section')

% checking chord lengths (should be equal)
initial_Length = x_Upper(1)-x_Lower(1);
transformed_Length = sqrt( (x_Upper_Rotated(1)-x_Lower_Rotated(1))^2+ (z_Rotated_Upper(1)- z_Rotated_Lower(1))^2 );

%% checking first profile
% checking first profile and curve matching in detail (top to bottom, same slope)
clf
plot([x_Upper_Rotated(first_y),x_Lower_Rotated(first_y)], [z_Rotated_Upper(first_y),z_Rotated_Lower(first_y)],'r');
hold on 
plot(first_Coordinates(:,1),first_Coordinates(:,3),'r');
scatter(first_Coordinates(:,1),first_Coordinates(:,3),'r');




%% Detialed transformation Steps


%more Detialed Testing (step by step of airfoil transformation) 
clf
chord_Length = x_Upper(last_y)-x_Lower(last_y);
transformed_Foil= fliped_Foil .* chord_Length;
plot(transformed_Foil(:,1),transformed_Foil(:,2),'r');
axis equal
hold on 
        
% translating to match initial guidance curve
transformed_Foil(:,1) = transformed_Foil(:,1) + x_Upper(last_y);
plot(transformed_Foil(:,1),transformed_Foil(:,2),'b');
axis equal

%now moving to quarter chord
transformed_Foil(:,1) = transformed_Foil(:,1) - quarterChord(last_y); 
plot(transformed_Foil(:,1),transformed_Foil(:,2),'g');
scatter(transformed_Foil(:,1),transformed_Foil(:,2),'g');
plot([-2,2],[0,0],'g'); 
axis equal


% now rotated
rotation = [cosd(twistAngle(last_y)), -sind(twistAngle(last_y)); sind(twistAngle(last_y)), cosd(twistAngle(last_y))];
transformed_Foil = transpose( rotation * transformed_Foil.'); 
plot(transformed_Foil(:,1),transformed_Foil(:,2),'y');
axis equal

 % now move back the "origin"
transformed_Foil(:,1) = transformed_Foil(:,1)- skewDistance(last_y)+ quarterChord(last_y);
plot(transformed_Foil(:,1),transformed_Foil(:,2),'k');
       
        
plot([x_Upper_Rotated(last_y),x_Lower_Rotated(last_y)], [z_Rotated_Upper(last_y),z_Rotated_Lower(last_y)],'b');
scatter(last_Coordinates(:,1),last_Coordinates(:,3),'b');
legend('stretched','matching Guide Curve','quarter Chord','quarterChordScatter','rotation','skewAndReorigin');

%% Testing Last Profile
% checking last profile and curve matching in detail (top to bottom, same slope)
clf
plot([x_Upper_Rotated(last_y),x_Lower_Rotated(last_y)], [z_Rotated_Upper(last_y),z_Rotated_Lower(last_y)],'b');
hold on 
plot(last_Coordinates(:,1),last_Coordinates(:,3),'b');
scatter(last_Coordinates(:,1),last_Coordinates(:,3),'b');
plot(transformed_Foil(:,1),transformed_Foil(:,2),'k');
 hold on
 
% also plotting initial maximum (x location) 
index = find(transformed_Foil(:,1) == max(transformed_Foil(:,1)));
scatter(transformed_Foil(index,1),transformed_Foil(index,2),'r'); 
 
 

  %% Trying to figure out why the point isnt on max
  
  % before transformation, making sure that the max point is on origin
  index = find(transformed_Foil(:,1) == max(transformed_Foil(:,1)));
  transformed_Foil(index,2) == 0
  
  % checking max point after rotation
  transformed_Foil(index,1) ==  x_Upper_Rotated(last_y)
  transformed_Foil(index,2) == z_Rotated_Upper(last_y)
    
    
    
 
  