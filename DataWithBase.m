function [combined_V , combined_F ,TR] = DataWithBase(Raw_data ,smoothing_factor)

% Reading the Data, smooth it, scale it and generate the faces and vertices

data = Raw_data

%turn the data table into array
data = table2array(data) ; 

%Extracting relative X and Y from Z matrix
[x, y] = meshgrid(1:size(data, 2), 1:size(data, 1));

%Lowering the number of data by decimation
decimationFactor = 2;
y_decimated = y(1:decimationFactor:end, :);
x_decimated = x(1:decimationFactor:end, :);
z_decimated = data(1:decimationFactor:end, :);

% Apply Gaussian smoothing
z_smoothed = imgaussfilt(z_decimated, smoothing_factor);

%Scaling
%This process have 2 parts, first defining an scale factor and then scaling
%scale factors for X and Y are defined in away that the final X and Y size would be 150 and 100
%the Z height is the mean of the X and Y
%X
x_ScaleFactor = 150 / max (x_decimated(:));
x_scaled = x_decimated * x_ScaleFactor;
%Y
y_ScaleFactor = 100 / max (y_decimated(:));
y_scaled = y_decimated * y_ScaleFactor;
%Z
z_ScaleFactor = ((100 + 150)/2) / max(z_smoothed(:)) ; %defining a scale factor to multiply all the Z and height to it (this transform the highest point of the graph to be mean of the X and Y size)%This make the printing easily 
z_scaled = z_smoothed * z_ScaleFactor ; %scale all the Z data 

% Flatten for triangulation
vertices = [x_scaled(:), y_scaled(:), z_scaled(:)]; %make vertices out of the data points(make a new matrix out of 3 matrix(each row is a 3D point
valid_indices = all(isfinite(vertices), 2); %remove unidentified values
vertices = vertices(valid_indices, :); %keeping the rows of the matrix with valid values

% generates a Delaunay triangulation using the X and Y coordinates from the vertices matrix
faces = delaunay(vertices(:, 1), vertices(:, 2));

% Create a solid volume from positive data
[F_solid, V_solid] = surf2solid(faces, vertices, 'elevation', -6); %2 is theheight of solid bottom(can be changed based on the data)
 



% Create Solid Cube

%creating a solid cube which the x dimension is 130 and y simension is 180
%the height of the cube is 5
vertices = [0 0 0;
    130 0 0;
    130 180 0;
    0 180 0;
    0 0 5;
    130 0 5;
    130 180 5;
    0 180 5];

%2D triangulation using x and y points
faces = delaunay(vertices(:, 1), vertices(:, 2));

%converting the surface into solid using surf2solid command
%put the elevation 10 and generating solid faces and solid vertices
[Solid_faces, Solid_vertices] = surf2solid(faces, vertices, 'elevation', 10);


% Add Custom Object on Top of Cube
% Position custom object centered on top of cube (e.g., 75, 100, 20)
custom_translation = [115, 15, 12]; 
rotation_xy = [cosd(90) -sind(90) 0; sind(90) cosd(90) 0; 0 0 1];
v_solid_rotated = (rotation_xy * V_solid')';
%moving 3D object
v_solid_translated = bsxfun(@plus, v_solid_rotated, custom_translation);
f_solid_shifted = F_solid + size(Solid_vertices, 1)


% Combine All and Export
combined_V = [Solid_vertices;v_solid_translated];
combined_F = [Solid_faces; f_solid_shifted];
% 
TR = triangulation(combined_F, combined_F);


% stlwrite(TR, 'cube_with_text_number_and_custom_object.stl');
% disp('✅ Exported as cube_with_text_number_and_custom_object.stl');

% Preview
% figure;
% patch('Faces', combined_F, 'Vertices', combined_V, ...
%       'FaceColor', [0.8 0.3 1], 'EdgeColor', 'none');
% axis equal;
% camlight;
% lighting gouraud;
% title('Cube with 3D Text, Numbers, and Custom Object');
% xlabel('X'); ylabel('Y'); zlabel('Z');
end %Data with base