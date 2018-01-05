%% Setting Up 
% Importing Master Data
filename = 'masterNormalized.txt' 
delimiterIn = '\t';
normalized = importdata(filename,delimiterIn);
cd 'P:\University\DuringUni\SubC\Fluids\props\RoughProps4\AirFoils\betterProfiles'


% Parameters:
profileNum = 10;
r_start = 0.03; %m
r_finish = 0.10; %m 


%{
% Parameters that needs to be changed
outputFileName = 'Section_0mm.txt'
chordL = 60* 1/100 ; % mm 
angle = 31.7/180 * pi;
zValue = 100; 
%}


%% Normalized and Translated Matrix
% testing code for transforming the original into a new and rotated cross
% section: 
    % 1: move the normalized airfoil so that the origin is now at the
    % center of lift
    normalized(:,1) = normalized(:,1)-1/4;
    plot(normalized(:,1), normalized(:,2))
    
    for n= 1: profileNum
        
        % rotation transformation
        z = ones(length(normalized(:,1)),1) .* n; 
        out = horzcat(rotate2D((normalized * 5*n).', 90).', z);
        
        % writing output
         dlmwrite(strcat('Section_',num2str(n),'.txt'),out,'delimiter','\t','precision',4)
    end
        
    
%% Functions 
% Rotation Matrix 
function rotatedMatrix = rotate2D(input, angle)
rotation = [cos(-angle), -sin(-angle); sin(-angle), cos(-angle)];
rotatedMatrix = rotation * input;
end

   