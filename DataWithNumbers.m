function[combined_V, combined_F , TR] = DataWithNumbers( Raw_data ,smoothing_factor , dimension_num)

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
[F_solid, V_solid] = surf2solid(faces, vertices, 'elevation', -6); %2 is theheight of solid bottom(can be changed based on the data)

% Create Solid Cube

%if the user put data for the file name od dimention values then goes to
%the

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

   % Add Multiple Numbers
%arrayfun : will do the function for each member of the input
%num2str : convert the numbers into characters
%'UniformOutput', false : we want the whole numbers to be in one cell not
%each store in a seperate cell so this is true when they all store in
%seperatedly cells
number_labels = arrayfun(@num2str, dimension_num, 'UniformOutput', false);
length(dimension_num(1))
number_positions = [
    102, 13, 10; %do not need adjusting
    102, 130, 10; 
    100,27, 10;
    0, 27, 10 %do not need adjusting
];

%adjust the position depending on the figure
if length(number_labels{2}) == 3
    number_positions(2 , :) = [102 , 140 , 10];
elseif length(number_labels{2}) == (1 || 2)
    number_positions(2 , :) = [102 , 150 , 10];
end
if length(number_labels{3}) == 2
    number_positions(3 , :) = [90 , 27 , 10];
elseif length(number_labels{3}) == 1
    number_positions(3 , :) = [90 , 27 , 10];
elseif length(number_labels{3}) == 3
    number_positions(3 , :) = [80 , 27 , 10];
end

%Rotates the object 270° around the Z-axis in XY plane
number_rotations = {eye(3), eye(3), ...
    [cosd(270) -sind(270) 0; sind(270) cosd(270) 0; 0 0 1], ...
    [cosd(270) -sind(270) 0; sind(270) cosd(270) 0; 0 0 1]};

num_vertices = [];
num_faces = [];
v_offset = 0;
fontSize = 15;
depth = 5;

for n = 1:length(number_labels)
    img = insertText(zeros(50,100), [0 0], number_labels{n}, ...
        'FontSize', fontSize,'Font', 'Arial Black', ...   % Use a naturally bold font
    'BoxOpacity', 0, ...
    'TextColor', 'white');
    BW = imbinarize(rgb2gray(img));
    BW = imrotate(BW, 270);
    BW = fliplr(BW);
    B = bwboundaries(BW, 'noholes');

    for k = 1:length(B)
        boundary = B{k}; % downsample
        x = boundary(:,2);
        y = boundary(:,1);
        z0 = zeros(size(x));
        z1 = depth * ones(size(x));
        coords0 = number_rotations{n} * [x'; y'; z0'];
        coords1 = number_rotations{n} * [x'; y'; z1'];
        minZ = min(coords0(3,:));
        coords0(3,:) = coords0(3,:) - minZ;
        coords1(3,:) = coords1(3,:) - minZ;
        coords0 = coords0 + number_positions(n, :)';
        coords1 = coords1 + number_positions(n, :)';
        numPts = size(coords0, 2);
        v0 = coords0';
        v1 = coords1';
        num_vertices = [num_vertices; v0; v1];
        for i = 2:numPts-1
            num_faces(end+1,:) = [v_offset+1, v_offset+i, v_offset+i+1];
            num_faces(end+1,:) = [v_offset+numPts+1, v_offset+numPts+i, v_offset+numPts+i+1];
        end
        for i = 1:numPts-1
            v1i = v_offset + i;
            v2i = v_offset + i + 1;
            v3i = v_offset + numPts + i + 1;
            v4i = v_offset + numPts + i;
            num_faces(end+1,:) = [v1i v2i v3i];
            num_faces(end+1,:) = [v1i v3i v4i];
        end
        v_offset = v_offset + 2 * numPts;
    end
end

%Add Custom Object on Top of Cube
% Assume user provides: f_solid, v_solid
% Position custom object centered on top of cube (e.g., 75, 100, 20)
custom_translation = [110, 20, 12]; 
rotation_xy = [cosd(90) -sind(90) 0; sind(90) cosd(90) 0; 0 0 1];
v_solid_rotated = (rotation_xy * V_solid')';
%moving 3D object
v_solid_translated = bsxfun(@plus, v_solid_rotated, custom_translation);
f_solid_shifted = F_solid + size(Solid_vertices, 1) + size(num_vertices, 1);

% Combine All and Export
num_faces_shifted = num_faces + size(Solid_vertices, 1);
combined_V = [Solid_vertices; num_vertices; v_solid_translated];
combined_F = [Solid_faces; num_faces_shifted; f_solid_shifted];
% 
TR = triangulation(combined_F, combined_V);
% stlwrite(TR , 'Cube with Numbers.stl');
% 
% figure;
% patch('Faces', combined_F, 'Vertices', combined_V, ...
%       'FaceColor', [0.8 0.3 1], 'EdgeColor', 'none');
% axis equal;
% camlight;
% lighting gouraud;
% title('Cube with Numbers');
% xlabel('X'); ylabel('Y'); zlabel('Z');
end
