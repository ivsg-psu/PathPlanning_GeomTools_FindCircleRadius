function [xc,yc,radii] = fcn_circleCenterFromThreePoints(x,y,varargin)
% FCN_circleCenterFromThreePoints calculates the center of a circle from
% three points given as vectors in x and y
% Format: 
% [xc,yc,radii] = fcn_circleCenterFromThreePoints(x,y,varargin)
%
% INPUTS:
%      x: a Nx1 vector where N is at least 3. If N = 1, a circle will be
%      fit between these threee points, if N = 4 or more, then one circle
%      will be fit to the first three points, another cicle to the next
%      three points, etc.
%
%      y: same dimension as x, but representing the y-coordinate of the
%      points
%
%      Optional input (3): a figure number if plot results are desired
%
% OUTPUTS:
%      xc: the x-coordinate of the centers of the circles, as an [(N-2)x1]
%      vector
%      yc: the y-coordinate of the centers of the circles, as an [(N-2)x1]
%      vector
%      radii: the radius of each the circles, as an [(N-2)x1]
%      vector
%
% Examples:
%      
%      % BASIC example
%      x = [0; 1; 0.5];
%      y = [0; 4; -1];
%      [xc,yc,radii] = fcn_circleCenterFromThreePoints(x,y,1)
% 
%      % ADVANCED example that uses vectors of x and y
%      x = [0; 1; 0.5; -1];
%      y = [0; 4; -1; 4];
%      [xc,yc,radii] = fcn_circleCenterFromThreePoints(x,y,1)
%
%      % ADVANCED example that lets user select N points 
%      figure(1); clf; grid on; axis equal;
%      [x,y] = ginput; % Get arbitrary N points until user hits return
%      [xc,yc,radii] = fcn_circleCenterFromThreePoints(x,y,1)     
%
%
% This function was written on 2020_03_20 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2020_03_20 - wrote the code
% 2020_05_22 - added more comments, particularly to explain inputs more
% clearly
%

do_debug = 1; % Flag to plot the results for debugging

%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Are there the right number of inputs?
if nargin < 2 || nargin > 3
    error('Incorrect number of input arguments')
end

% Does user want to show the plots?
if 3 == nargin
    fig_num = varargin{1};
    figure(fig_num);
    do_debug = 1;
else
    if do_debug
        fig = figure; %#ok<NASGU> % create new figure with next default index
    end
end

% Are the input vectors the right shape?
if length(x(:,1))<3
    error('x vector must be at least a 3x1');
end
if length(y(:,1))<3
    error('y vector must be at least a 3x1');
end

%% Solve for the circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% The code below is set up to be vectorized if there are more than one
% solution. Since the code is quite different looking for each, they are
% separated out. However, it may be that the N-solution case works for N is
% equal to 1. This was not tested.

% Do some pre-calculations
num_solutions = length(x(:,1))-2; % This is the number of solutions to expect
r_squared = (x.^2 + y.^2); % These are the radii of points from origin
diff_x = diff(x);
diff_y = diff(y);
diff_rsquared = diff(r_squared);

if 1 == num_solutions % Expecting just one solution. No need for big A, b matrices    
    % solve for the center point    
    A = [diff_x diff_y];
    b = 1/2*diff_rsquared;
    
else % Simultaneous solutions to be calculated - create big A and b matrices    
    % Construct the A-matrix and b matrix that will create the regressor.
    % Start by filling A and b matrices up with zeros (see notes for
    % explanation of iputs)
    A = zeros(2*num_solutions,2*num_solutions);
    b = zeros(2*num_solutions,1);
    
    % Fill in the non-zero portions of the matrix, which will be 1 per each
    % of the N solutions
    for i_solution = 1:num_solutions
        A(1+2*(i_solution-1):2+2*(i_solution-1),1+2*(i_solution-1)) = ...
            diff_x(i_solution:i_solution+1);
        A(1+2*(i_solution-1):2+2*(i_solution-1),2+2*(i_solution-1)) = ...
            diff_y(i_solution:i_solution+1);        
        b(1+2*(i_solution-1):2+2*(i_solution-1),1) = ...
            1/2*diff_rsquared(i_solution:i_solution+1);
    end
end

% Solve the center points
centers = A\b;
centers = reshape(centers,2,length(centers(:,1))/2);
xc = centers(1,:)';
yc = centers(2,:)';


% NOTE: the following line is the slowest in the code. It can be sped
% up if we do not take the square root
radii = ((xc - x(1:num_solutions,1)).^2 + ...
    (yc - y(1:num_solutions,1)).^2).^0.5;

%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_debug
    hold on % allow multiple plot calls
    plot(x,y,'ro');  % Plot all the input points
    plot(xc,yc,'g+'); % Plot all the circle centers

    axis equal;
    grid on; grid minor;

    % plot all the circle fits
    angles = 0:0.01:2*pi;
    for i_fit = 1:length(xc)       
      
        x_circle = xc(i_fit,1) + radii(i_fit) * cos(angles);
        y_circle = yc(i_fit,1) + radii(i_fit) * sin(angles);
        plot(x_circle,y_circle,'b-');
    end
end



