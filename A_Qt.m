function A = A_Qt(Q, t)
    N = length(Q(1,1,:));

    A = zeros(3);
    for i = 1:N
        A = A + Q(:,:,i) * t;
    end
end