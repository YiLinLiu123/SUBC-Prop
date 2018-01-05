%% Setting Up 
% Importing Master Data
filename = 'masterCoordinates.txt';
delimiterIn = '\t';
original = importdata(filename,delimiterIn);
cd 'P:\University\DuringUni\SubC\Fluids\props\RoughProps4\AirFoils\profiles'


% Parameters that needs to be changed
outputFileName = 'Section_0mm.txt'
chordL = 60* 1/100 ; % mm 
angle = 31.7/180 * pi;
zValue = 100; 


% plotting 
plot(original(:,1), original(:,2))



%% Normalized and Translated Matrix
% testing code for transforming the original into a new and rotated cross
% section: 
normalized  = original ./60; % normalized down to chord of 1mm

    % 1: move the normalized airfoil so that the origin is now at the
    % center of lift
    normalized(:,1) = normalized(:,1)-1/4;
    
    % 2: Scale the normalized airfoil and then rotate to desired angle
    z = ones(length(normalized(:,1)),1) .* zValue; 
    out = horzcat(rotate2D((normalized * chordL).', angle).', z);
    
    % 3: save to the output file:
    dlmwrite(outputFileName,out,'delimiter','\t','precision',4)
    
    
    plot(out(:,1),out(:,2))   
    pwd
    
%% Functions 
% Rotation Matrix 
function rotatedMatrix = rotate2D(input, angle)
rotation = [cos(-angle), -sin(-angle); sin(-angle), cos(-angle)];
rotatedMatrix = rotation * input;
end

   