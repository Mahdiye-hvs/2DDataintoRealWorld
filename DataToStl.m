%%GC×GC subfunction

function [F_solid, V_solid , tri] = DataToStl( Raw_data ,smoothing_factor )
%This functions turn a raw data into Stl file, ready for print

%inputs:
%Raw_data : one column data, in Csv ot txt format
%smoothing_factor

%outputs:
%3D plot of the solid object
%Stl file saved as GC×GC.stl

% Read the file
data = Raw_data;

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
[F_solid, V_solid] = surf2solid(faces, vertices, 'elevation', -5); %2 is theheight of solid bottom(can be changed based on the data)

% Save STL file for positive peaks
tri = triangulation(F_solid, V_solid); %This command creates a triangulation object in MATLAB, which provides geometric and topological information about the mesh.

end




 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 