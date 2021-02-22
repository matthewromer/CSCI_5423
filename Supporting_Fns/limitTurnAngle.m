function [limitedDir] = limitTurnAngle(currentDir,targetDir,maxAngle)
%limitTurnAngle updates the target direction based on a maximum turn angle
%   Inputs:  Current direction, initial target direction, maximum turn
%            angle (radians)
%   Outputs: Updated target direction with turn angle limited 
%
%   Notes:   Method of computing angle direction comes from a stack
%            overflow question: https://stackoverflow.com/q/27681893
%   

%Compute angle between current angle and origin and target angle and origin
currAng = atan2(currentDir(2),currentDir(1));
targAng = atan2(targetDir(2),targetDir(1));

%Compute difference in angles 
deltaAng = acos(dot(currentDir,targetDir)/norm(currentDir)/norm(targetDir));

%If angle is greater than the maximum, change the target direction to one
%that is no more than the maximum angle different from the initial
%direction
if deltaAng > maxAngle
    
    %Determine whether to turn CW or CCW
    q = sin(currAng-targAng);
    if q > 0
        dir = -1;
    else
        dir = 1;
    end
        
    %Form rotation matrix and turn target direction 
    theta = maxAngle*dir;
    mat = [cos(theta) -sin(theta);
           sin(theta)  cos(theta)];
         
    limitedDir = mat*currentDir;
else
    limitedDir = targetDir;
end

end