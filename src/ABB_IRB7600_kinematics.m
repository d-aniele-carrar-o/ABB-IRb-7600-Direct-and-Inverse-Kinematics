clc;
close all;

%%
scaleFactor = 1;

% if true, compute AND plot the manipulator configuration
draw = true;

% axis limits for IRB 7600-400/2.55, 404 mm for LeanID (not used yet)
limits = [deg2rad(180), deg2rad(-180);  % q1
          deg2rad(85),  deg2rad(-60);   % q2
          deg2rad(60),  deg2rad(-180);  % q3
          deg2rad(300), deg2rad(-300);  % q4
          deg2rad(100), deg2rad(-100);  % q5
          deg2rad(220), deg2rad(-220);  % q6
          ];

% D-H parameters:
% T   :    0->1  1->2   2->3    3->4    4->5    5->6
% i   = | 0 |  1  |   2  |   3   |   4   |   5   |   6   |
alpha = [  0, pi/2,     0,   pi/2,  -pi/2,   pi/2,     0 ];
a     = [  0, 0.41, 1.075,  0.165,      0,      0,     0 ];
d     = [  0, 0.78,     0,      0,  1.056,      0,  0.25 ];

% joint angles configuration
q = [0.33, 2.476, -1.189, 2.127, 0.563, -2.138]

% compute direct kinematics for defined joint angles q
T06 = direct_kinematics( q, alpha, a, d, limits )

% compute 8 (at most) inverse kinematics solutions for given Frame_6 pose
H = inverse_kinematics( T06, alpha, a, d )

% for every inverse kinematics solution, compute direct kinematics and draw the manipulator pose
for i=1:8
    i
    q_i = H(i,:);
    if anynan(q_i)
        fprintf("configuration impossible, moving to next one.")
        continue
    end
    if draw
        if i==1
            lim = 1;
            limS = scaleFactor*lim;
            alfa = 340;
            beta = 140;
            T06 = direct_kinematics( q_i, alpha, a, d, limits )
            % axs  = axes( 'XLim', [-limS/2 limS], 'YLim', [-limS/2 limS/2], 'ZLim', [-0.1 limS] );
            l = max([T06(1,4), T06(2,4), T06(3,4)]);
            axs  = axes( 'XLim', [-1.5*l, 1.5*l], 'YLim', [-1.5*l, 1.5*l], 'ZLim', [-0.1, 1.5*l] );
            view( alfa, beta ); grid on;
            xlabel(['X x ', num2str(scaleFactor)], 'FontSize', 12);
            ylabel(['Y x ', num2str(scaleFactor)], 'FontSize', 12);
            zlabel(['Z x ', num2str(scaleFactor)], 'FontSize', 12);
            handles(1) = axs;
            [T06, handlesR, frame_names] = direct_kinematics_draw( q_i, alpha, a, d, limits, handles, true, scaleFactor );
        else
            [T06, handlesR, frame_names] = direct_kinematics_draw( q_i, alpha, a, d, limits, handlesR, false, scaleFactor );
            T06
        end
        pause()
        
        % clear frame names labels
        if i < 8
            delete(frame_names)
        end
    else
        T06 = direct_kinematics( q_i, alpha, a, d, limits )
    end
end


