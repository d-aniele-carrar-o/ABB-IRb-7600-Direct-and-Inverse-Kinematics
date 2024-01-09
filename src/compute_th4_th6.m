% solve non-linear system of equations
function [th4, th6] = compute_th4_th6( T36_31, T36_32, th5 )
    % | gamma = th4 + th6
    % | sin(gamma) = sin(th4)*cos(th6) + cos(th4)*sin(th6)
    % | cos(gamma) = cos(th4)*cos(th6) - sin(th4)*sin(th6)
    
    F = @(th4, th6) [ cos(th4)*sin(th6) + cos(th5)*sin(th4)*cos(th6);
                      cos(th4)*cos(th6) - cos(th5)*sin(th4)*sin(th6) ];
    C = [ T36_31; T36_32 ];

    f = @(u) F( u(1), u(2) ) - C;

    u0 = [0.5; 0.5];

    u = fsolve( f, u0 );

    r = f(u);

    th4 = r(1);
    th6 = r(2);
end
