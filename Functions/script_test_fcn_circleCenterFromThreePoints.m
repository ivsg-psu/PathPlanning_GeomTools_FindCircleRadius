% script_test_fcn_circleCenterFromThreePoints - this is is a script written
% to test the function: fcn_circleCenterFromThreePoints.m.
%
% Revision history:
% 2020_03_20 - started writing the function and script
% 2020_05_22 - added more comments

x = [0; 0.5; 1; 4; 6; 7; 9; 11; 15];
y = [0; 4;  -1;-3; 2; -1;3;  3; -0.5];
fcn_circleCenterFromThreePoints(x,y,1);
plot(x,y,'r-');
hold on
figure(1); clf;
for i=1:length(x)-2
    fcn_circleCenterFromThreePoints(x(i:i+2),y(i:i+2),1);
    plot(x,y,'g-');
    %pause;
end


% x = [0; 1; 0.5; 5];
% y = [0; 4; -1; 6];
% fcn_circleCenterFromThreePoints(x,y,1);


% The following tests the N input form:
button = 1;
while sum(button) <=1   % read ginputs until a mouse right-button occurs   
    % Get a new point and redo plot
    [x(end+1),y(end+1),button] = ginput(1); %#ok<SAGROW>
    fcn_circleCenterFromThreePoints(x,y,1);     
end


% % The following tests the 3 input form:
% button = 1;
% while sum(button) <=1   % read ginputs until a mouse right-button occurs
%     % Shift points down to prep for next input
%     x(1:end-1) = x(2:end);
%     y(1:end-1) = y(2:end);
%     
%     % Get a new point and redo plot
%     [x(end),y(end),button] = ginput(1);
%     fcn_circleCenterFromThreePoints(x,y,1);     
% end