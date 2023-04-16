function D = D_Qt(Q, t)
    N = length(Q(1,1,:));
    
    h = N*t;

    D = zeros(3,3);
    for i = 1:N
        zk = -h/2 + i*t;
        zk_1 = -h/2 + (i-1)*t;

        D = D + Q(:,:,i)/3 * (zk^3 - zk_1^3);
    end
end