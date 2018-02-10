%% Parameters you input

clear
% files will be saved at which directory specified here. Please modify to
% desired folder path
cd 'P:\University\DuringUni\SubC\Fluids\props\RoughProps4\AirFoils\ReorganizedCode'

% output file format
out_File = 1; %1 for txt (tab delimited), 0 for xlsx
leftOrRight = -1; % 1 for left handed, -1 for right handed

% airfoil Related:
airfoil_Filename = 'masterNormalized.txt'; % assuming the airfoil is normalized
delimiterIn = '\t'; % delimiter type of the airfoil file name

% initial guide curve related: (the ellipses that was given to you)
startY = 0; % (cm)
endY = 10; %ending Y (cm)
points = 300+1; %number of data points you want, add one to account for the fact that the startin point is included
y = linspace(startY,endY,points); %% radius direction, generating the list of y-coordinates 
x_Upper  = ((13^2 - (2.*y-7).^2)./36).^(1/2);
x_Lower = ((13^2 - (2.*y-7).^2)./100).^(1/2)* -1;


%twist angle function: the function from stephen's excel sheets
a_4 = 178767; % fourth power coefficient
a_3 = 87784; % 3rd power coefficient
a_2 = 17140; %2nd power coefficient
a_1 = 1721.7;
a_0 =100.11;

% airfoil profiles
% Generating the first profile (at the "origin")
first_y = 1/5*(points-1); %point number please, location of first profile (mm along radius)
last_y = points-1; % point number please


%% Behind the scene 
% importing airFoil
normalizedFoil = importdata(airfoil_Filename,delimiterIn);

% 1/4 chord locations:
quarterChord = x_Lower + 3/4*(x_Upper-x_Lower) +0.1*(max(y)-y).^(1.1);
% [zeros(1,extra_Mod), .01* linspace(extra_Mod, points,(points-extra_Mod))];
quarterChord_Orig = x_Lower + 3/4*(x_Upper-x_Lower);
% setting up the airfoil:
fliped_Foil= normalizedFoil;
fliped_Foil(:,1) = -fliped_Foil(:,1); % flipping the airfoil about y-axis
fliped_Foil(:,2) = fliped_Foil(:,2); % flipping the airfoil about y-axis

%% Testing for Airfoil
clf
 
% plotting the airfoil
title('NACA-2411 Airfoil')
plot(normalizedFoil(:,1),normalizedFoil(:,2), 'r');
hold on
axis equal

plot(fliped_Foil(:,1), fliped_Foil(:,2),'b'); 
scatter(fliped_Foil(:,1), fliped_Foil(:,2),'b');
plot([-2,2],[0,0]); 
hold off


%% Test initial guide Curves
% plotting the initial guidance curves
clf
title('Flat Propellor Curves')
plot(y,x_Upper,'r');
hold on 
plot(y,x_Lower,'b');
plot(y,quarterChord,'y');
plot(y,quarterChord_Orig);
xlabel('y-coordinates (mm)') % x-axis label
ylabel('x-coordinates (mm)')
axis equal


%% Transformations

twistAngle = a_4*(y./100).^4-a_3*(y./100).^3+ a_2*(y./100).^2-a_1*(y./100)+a_0;
skewAngle = 0.1 * (y-startY).^(1.5); %function for skew angle along radius (degrees)
skewDistance = tand(skewAngle) .* (y-startY); % calculating the skew distance along y-axis (mm)

%Hotfix, subtract a bit of x


% place holders for curve transformations:
x_Upper_Rotated = zeros(1,points);
x_Lower_Rotated = zeros(1,points);
z_Rotated_Upper = zeros(1,points); 
z_Rotated_Lower = zeros(1,points);



% master for loop tranformations
for n= 1: points
      % the manuipulation is done for each x-z point
    target = [x_Upper(n),x_Lower(n);0,0]; 
    
    % moving origin to 1/4 chord length (z still 0) 
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
        transformed_Foil(:,1) = transformed_Foil(:,1);
        %take apart some of the z-component
        transformed_Foil(:,2) = leftOrRight* transformed_Foil(:,2);
        
        % storing information
        if(n == first_y)
           first_Profile = transformed_Foil;  
        elseif (n== last_y)
           last_Profile = transformed_Foil; 
        end
        
        
    end 
    
    
    
end 

%% writing out the curve files (txt files)

z_Rotated_Upper =leftOrRight*z_Rotated_Upper;
z_Rotated_Lower = leftOrRight* z_Rotated_Lower; 

% reorganizing the matrices. 
twisted_Upper = vertcat(x_Upper_Rotated,y-startY, z_Rotated_Upper); 
twisted_Lower = vertcat(x_Lower_Rotated,y-startY,z_Rotated_Lower);

% generating the corresponding "y" components
first_Profile_y = zeros(1, length(first_Profile))+y(first_y)-startY; 
last_Profile_y =zeros(1, length(last_Profile))+y(last_y)-startY;

%formatting profiles
 first_Coordinates = horzcat(first_Profile(:,1), transpose(first_Profile_y), first_Profile(:,2));
 last_Coordinates = horzcat(last_Profile(:,1), transpose(last_Profile_y), last_Profile(:,2));

 
 if(out_File == 0)
    if(leftOrRight == 1)
        xlswrite('twisted_Upper_Left.xlsx',transpose(twisted_Upper));
        xlswrite('twisted_Lower_Left.xlsx',transpose(twisted_Lower)); 
        xlswrite('first_Profile_Left.xlsx',first_Coordinates);
        xlswrite('last_Profile_Left.xlsx',last_Coordinates);
     elseif(leftOrRight == -1)
        xlswrite('twisted_Upper_Right.xlsx',transpose(twisted_Upper));
        xlswrite('twisted_Lower_Right.xlsx',transpose(twisted_Lower)); 
        xlswrite('first_Profile_Right.xlsx',first_Coordinates);
        xlswrite('last_Profile_Right.xlsx',last_Coordinates);
    end
 elseif (out_File== 1)
     if(leftOrRight == 1)
         dlmwrite('twisted_Upper_Left.txt',transpose(twisted_Upper),'delimiter','\t','precision',5);
         dlmwrite('twisted_Lower_Left.txt',transpose(twisted_Lower),'delimiter','\t','precision',5);
         dlmwrite('first_Profile_Left.txt',first_Coordinates, 'delimiter','\t','precision',5);
         dlmwrite('last_Profile_Left.txt',last_Coordinates, 'delimiter','\t','precision',5);
     elseif(leftOrRight == -1)
         dlmwrite('twisted_Upper_Right.txt',transpose(twisted_Upper),'delimiter','\t','precision',5);
         dlmwrite('twisted_Lower_Right.txt',transpose(twisted_Lower),'delimiter','\t','precision',5);
         dlmwrite('first_Profile_Right.txt',first_Coordinates, 'delimiter','\t','precision',5);
         dlmwrite('last_Profile_Right.txt',last_Coordinates, 'delimiter','\t','precision',5);
     end 
 end

 
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
