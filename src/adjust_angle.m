function a=adjust_angle(a)
    %{
        Adjust angle to accomidate minimal angle measurement methods
        This reports a minimal angle that is not specific to the vertical
        or horizontal axis.
    %}

    if a > 45 && a <= 90
        a = 90 - a;
    elseif a > 90 && a <= 135
        a = a - 90;
    elseif a > 135 && a <= 180
        a = 180 - a;
    elseif a < -135 && a >= -180
        a = a * -1;
        a = 180 - a;
    elseif a < -90 && a >= -135
        a = a * -1;
        a = a - 90;
    elseif a < -45 && a >= -90
        a = a * -1;
        a = 90 - a;
    elseif a < 0 && a >= -45
        a = a * -1;
    end
end