function [x,y]=find_ver_hor(angle, h)
	%
	%	Translate angle and distance into proportion of horizonal and 
	%  vertical movement needed to create movement in the direction of the angle.
	%  Distance is hypotenuse of right triangle with ver (vertical rise)
	%  and hor (horizontal run) as sides. Given angle and hypotenuse the
	%  sides are found using the Pythagorean theorem. h = hypotenuse of triangle.
	%  Reference: https:%www.mathsisfun.com/algebra/trig-finding-side-right-triangle.html
	%
	angle = (angle/360)*2*pi; % convert from degrees to radians

	if (angle < (pi/2))
		x = sin(angle) * h;
		y = sqrt((h^2)-(x^2));
    elseif (angle >= (pi/2) && angle < pi)
		x = cos(angle-(pi/2)) * h;
		y = sqrt((h^2)-(x^2)) * -1;
    elseif (angle >= pi && angle < (pi*1.5))
		x = cos((pi*1.5)-angle) * h * -1;
		y = sqrt((h^2)-(x^2)) * -1;
    elseif (angle >= (pi*1.5) && angle <= (pi*2))
		x = cos(angle-(pi*1.5)) * h * -1;
		y = sqrt((h^2)-(x^2));
	end
end