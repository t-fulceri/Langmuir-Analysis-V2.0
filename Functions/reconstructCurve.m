function [ output_args ] = reconstructCurve( x_of_t, y_of_t)
%RECONSTRUCT_CURVE Reconstructs the y = y(x) curve starting from the x(t) and
%y(t) signals, taking x as the independent variable
%   Warning: x and y must be of the same length and must represent
%   timeseries of measurements taken simultaneously

[x_axis, unique_ind_array] = unique(x_of_t);
y_of_x = y_of_t(unique_ind_array);

for i = 1:length(unique_ind_array)
    unique_ind = unique_ind_array(i);
    %Find the set of indices of x_of_t which correspond to the same x value
    same_x_ind_array = find(x_of_t == x_of_t(unique_ind));
    y_values = y_of_t(same_x_ind_array);
    y_of_x(unique_ind) = mean(y_values);
end
y_of_x = y_of_x(unique_ind_array);

%Copy the curve and the uncertainty into the output structure
output_args.x = x_axis;
output_args.y_of_x = y_of_x;

end
