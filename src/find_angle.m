function a=find_angle(y2,y1,x2,x1)
    %a=atand(abs(y2-y1)/abs(x2-x1)); % angle in degrees
    a=atan2(y2 - y1, x2 - x1) * 180 / pi; % angle in degrees
end