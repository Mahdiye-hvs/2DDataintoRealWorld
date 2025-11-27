function [combined_V , combined_F ,TR] = DataWithBase(Raw_data ,smoothing_factor)

%inputs:
%Raw_data : one column data, in Csv ot txt format
%smoothing_factor

%outputs:
%combined_V: Faces generated from the data points(including the added height base)
%combined_F: Vertices generated from the data points(including the added height base)
%TR: Triangulated objects ready for exporting the STL format(including the added height base)


% Reading the data
data = Raw_data

% Turn the data table into array
data = table2array(data) ;

% Extracting relative X and Y from Z matrix
[x, y] = meshgrid(1:size(data, 2), 1:size(data, 1));

% Lowering the number of data by decimation
decimationFactor = 2;
y_decimated = y(1:decimationFactor:end, :);
x_decimated = x(1:decimationFactor:end, :);
z_decimated = data(1:decimationFactor:end, :);

% Apply Gaussian smoothing
z_smoothed = imgaussfilt(z_decimated, smoothing_factor);

% Scaling
% This process have 2 parts, first defining an scale factor and then scaling
% Scale factors for X and Y are defined in away that the final X and Y size would be 150 and 100
% The Z height is the mean of the X and Y
% X
x_ScaleFactor = 150 / max (x_decimated(:));
x_scaled = x_decimated * x_ScaleFactor;
% Y
y_ScaleFactor = 100 / max (y_decimated(:));
y_scaled = y_decimated * y_ScaleFactor;
% Z
z_ScaleFactor = ((100 + 150)/2) / max(z_smoothed(:)) ; %defining a scale factor to multiply all the Z and height to it (this transform the highest point of the graph to be mean of the X and Y size)%This make the printing easily
z_scaled = z_smoothed * z_ScaleFactor ; %scale all the Z data

% Flatten for triangulation
vertices = [x_scaled(:), y_scaled(:), z_scaled(:)]; %make vertices out of the data points(make a new matrix out of 3 matrix(each row is a 3D point
valid_indices = all(isfinite(vertices), 2); %remove unidentified values
vertices = vertices(valid_indices, :); %keeping the rows of the matrix with valid values

% Generates a Delaunay triangulation using the X and Y coordinates from the vertices matrix
faces = delaunay(vertices(:, 1), vertices(:, 2));

% Create a solid volume from positive data
[F_solid, V_solid] = surf2solid(faces, vertices, 'elevation', -6); %2 is theheight of solid bottom(can be changed based on the data)

% Create solid cube
% Creating a solid cube which the x dimension is 130 and y simension is 180
% The height of the cube is 5
vertices = [0 0 0;
    130 0 0;
    130 180 0;
    0 180 0;
    0 0 5;
    130 0 5;
    130 180 5;
    0 180 5];

% 2D triangulation using x and y points
faces = delaunay(vertices(:, 1), vertices(:, 2));

% Converting the surface into solid using surf2solid command
% Put the elevation to 10 and generating solid faces and solid vertices
[Solid_faces, Solid_vertices] = surf2solid(faces, vertices, 'elevation', 10);


% Add custom object on top of cube
% Position custom object centered on top of cube
custom_translation = [115, 15, 12];
rotation_xy = [cosd(90) -sind(90) 0; sind(90) cosd(90) 0; 0 0 1];
v_solid_rotated = (rotation_xy * V_solid')';
% Moving 3D object and matching the faces
v_solid_translated = bsxfun(@plus, v_solid_rotated, custom_translation);
f_solid_shifted = F_solid + size(Solid_vertices, 1);

% Combine all and export
combined_V = [Solid_vertices;v_solid_translated];
combined_F = [Solid_faces; f_solid_shifted];

% Build triangulation
TR = triangulation(combined_F, combined_F);

end %DataWithBase
