function sig = rotate_stress(sig, theta)
    theta= theta*2*pi/360;
    c = cos(theta);
    s = sin(theta);
            
    Ts = [c^2       s^2     2*c*s
          s^2       c^2     -2*c*s
          -c*s      c*s     c^2-s^2]; 

    sig = Ts * sig;
end