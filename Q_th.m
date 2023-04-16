% code from Giovanni
function Q = Q_th(E1, E2, nu12, G12, thetas)
    nu21 = nu12 / E1 * E2;
    
    Q = zeros(3, 3, length(thetas));
    for i = 1:length(thetas)
        th = thetas(i);
    
        Q12 =   [E1/(1-nu12*nu21)       nu12*E2/(1-nu12*nu21)   0
                nu12*E2/(1-nu12*nu21)   E2/(1-nu12*nu21)        0
                0                       0                       G12];
        
        c = cosd(th);
        s = sind(th);
          
        Ts = [c^2       s^2     2*c*s
              s^2       c^2     -2*c*s
              -c*s      c*s     c^2-s^2];    
        Te = [c^2       s^2     c*s
              s^2       c^2     -c*s
              -2*c*s    2*c*s   c^2-s^2];
        Q(:,:,i)  = Ts\Q12*Te;
    end
end