function [line_handle, poly_handle] =  plot_filled_ci(x, y, CI, colorTrip)
% function [line_handle, poly_handle] =  plot_filled_ci(x, y, CI, colorTrip)
%
% PURPOSE:
%   This function is an edit of JB Trimper's code "error_fill_plot2" that
%   creates an error fill plot using confidence intervals, so that the
%   filled error can vary + and - at each point.
%
% INPUT:
%           x = 1 x n x values
%           y = 1 x n y values
%           CI = 2 x n CI values
%   colorTrip = an RGB triplet specifying the color to use
%
% OUPUT:
%   line_handle = handle for the line that is plotted
%   poly_handle = handle for the shaded error polygon
%
% MM Donahue - edited from JB Trimper
% 03/2020 - edited from 10/20/14
% Colgin Lab - edit from Manns Lab

if length(x)~=length(y) || length(x)~=length(CI) || length(y)~=length(CI)
    error 'Length of inputs must be the same.';
end

if (size(x,2)==1); x=x'; end
if (size(y,2)==1); y=y'; end
if (size(CI,2)==2); CI=CI'; end

upper_error= CI(2,:);
lower_error= CI(1,:);

x_poly=[x, fliplr(x)];
hold on;
error_poly=[lower_error, fliplr(upper_error)];

poly_handle=fill(x_poly,error_poly, colorTrip);
set(poly_handle, 'edgecolor', colorTrip);
alpha(.3);

line_handle=plot(x,y,'color', colorTrip, 'LineWidth', 2);


end %fnctn

