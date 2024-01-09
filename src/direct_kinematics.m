function [T06] = direct_kinematics( q, alpha, a, d, limits )
    % homogeneous transformation matrix from frame [i-1] to frame [i]:
    % [i-1]T[i] = Rot(x[i-1],alpha(i-1)) * Trasl(x[i-1],a[i-1]) * Rot(z[i],theta[i]) * Trasl(z[i],d[i])

    % homogeneous transformation matrix frame 0 -> 1
    T01 = @(q1) rot_tralsX( alpha(1), a(1) ) * rot_traslZ( q1, d(2) );
    
    % homogeneous transformation matrix frame 1 -> 2
    T12 = @(q2) rot_tralsX( alpha(2), a(2) ) * rot_traslZ( q2, d(3) );
    
    % homogeneous transformation matrix frame 2 -> 3
    T23 = @(q3) rot_tralsX( alpha(3), a(3) ) * rot_traslZ( q3, d(4) );
    
    % homogeneous transformation matrix frame 3 -> 4
    T34 = @(q4) rot_tralsX( alpha(4), a(4) ) * rot_traslZ( q4, d(5) );

    % homogeneous transformation matrix frame 4 -> 5
    T45 = @(q5) rot_tralsX( alpha(5), a(5) ) * rot_traslZ( q5, d(6) );

    % homogeneous transformation matrix frame 5 -> 6
    T56 = @(q6) rot_tralsX( alpha(6), a(6) ) * rot_traslZ( q6, d(7) );

    % check joint limits
    % q = check_limits( q, limits );

    % compute matrices with actual joint variables values q=[q1, q2, q3, q4, q5, q6]
    T01m   = T01(q(1));
    T12m   = T12(q(2));
    T23m   = T23(q(3));
    T34m   = T34(q(4));
    T45m   = T45(q(5));
    T56m   = T56(q(6));
    
    % compute total transformation matrix frame 0 -> 6
    T06 = T01m*T12m*T23m*T34m*T45m*T56m;
end

