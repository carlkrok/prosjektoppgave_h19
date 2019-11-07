function rotMat = RotMatXYZEulerDeg(rollDeg, pitchDeg, yawDeg)

    Rx =    [ 1, 0, 0;
            0, cosd(rollDeg), -sind(rollDeg);
            0, sind(rollDeg), cosd(rollDeg)];
        
    Ry =    [ cosd(pitchDeg), 0, sind(pitchDeg);
            0, 1, 0;
            -sind(pitchDeg), 0, cosd(pitchDeg)];
        
    Rz =    [ cosd(yawDeg), -sind(yawDeg), 0;
            sind(yawDeg), cosd(yawDeg), 0;
            0, 0, 1];
        
    rotMat = Rz * Ry * Rx;
    
end

