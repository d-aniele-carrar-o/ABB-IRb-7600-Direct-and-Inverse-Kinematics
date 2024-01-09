% homogeneous transformation matrix from frame [i-1] to frame [i]
function [T] = transf_i_1_i( i, theta_i, alpha, a, d )
    T = rot_tralsX( alpha(i), a(i) ) * rot_traslZ( theta_i, d(i+1) );
end
