
function [combined_V, combined_F , TR] = DataWithDetails( Raw_data ,smoothing_factor , Data_name , dimension_num)
%This functions turn a raw data into Stl file, ready for print
%the final STL file will have a name of the file, and the numbers related
%to the dimensions

%inputs:
%Raw_data : one column data, in Csv ot txt format
%smoothing_factor
%Data_name
%dimension_num : will contains 4 number as an array which each is related
%to the starting and ending point of the x and y data(while the z data is
%about the detctor's intensity, 
% example: {120 , 240 , 350 , 780}
%120 : y starting point
%240 : y ending point
%350 : x starting point
%780 : x ending point

%outputs:
%3D plot of the solid object
%Stl file saved with the number and text on it
%% Reading the Data, smooth it, scale it and generate the faces and vertices

data =Raw_data;

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
 
%% Create Solid Cube

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

%% Create 3D Text Letter by Letter 

%get the input letters
letters = char(Data_name);
fontSize = 15; %font size
depth = 5; %the amount of extruding
spacing = 12; %the spacing between each character
startY = 0; %the starting point

%for making the letters as standing ones we need a matrix to multiply and
%turn the upright
R = [cos(pi/2) 0 sin(pi/2); 0 1 0; -sin(pi/2) 0 cos(pi/2)];

%make an empty cells for the words
%creates a cell which has 1 row and length * 2 column
%*2 = 2 vertices, one for back and one for front(like faces)
all_text_vertices = cell(1, length(letters)*2);
all_text_faces = cell(1, length(letters)*2);

%how many vertices have already been added 
vertex_offset = 0;

%Track where you are in your list of geometry pieces 
face_idx = 1;


for i = 1:length(letters)

    if letters(i) == ' ' %check if the character is blank
        startY = startY + spacing;%change the starting point

        %if the characters is space, then when its reached to this line it
        %will start over a new iteration in for loop
        continue 
    end
    
    %turn the text into numbers and matrix(its not binary yet)
    %generate a RGB one
    %generate a 50 × 50 gray scale which the insert text will return it
    %into RGB
    %starting point is (0 , 0)
    %adjusting the font and font size
    %turn the background opocity to zero
    %the text colors are white
    img = insertText(zeros(50,50), [0 0], letters(i), ...
         'FontSize', fontSize,'Font', 'Arial Black', ...   % Use a naturally bold font
    'BoxOpacity', 0, ...
    'TextColor', 'white');

    %make the RGB into gray + make the image binary
    BW = imbinarize(rgb2gray(img));
    %rotating 270 degree and flip it
    %this is essential for showing the letters
    BW = imrotate(BW, 270);
    BW = fliplr(BW);
    %extract the boundary coordinates of white shapes
    B= bwboundaries(BW, 'noholes'); 
    %placing a letters in 130 x, change the starting point for the next
    %letter , and Z in 0
    pos = [130; 10 + startY; 0];

    %looping through all detected boundaries in the binary image 
    for k = 1:length(B)
        % reduces the number of points
        boundary = B{k}; % downsample for speed
        
        %extracting x and y
        x = boundary(:,2);
        y = boundary(:,1);
        z0 = zeros(size(x)); %create zero z matrix
        z1 = depth * ones(size(x)); %create the matrix of z with depth
        %the bottom face
        %make a matrix of x,y,and z which each are in a row
        %multiply to a matrix which make it stands upright
        coords0 = R * [x'; y'; z0'];
        %the top face
        coords1 = R * [x'; y'; z1'];
        minZ = min(coords0(3,:));%find min z
        coords0(3,:) = coords0(3,:) - minZ; %shifting z to be zero
        coords1(3,:) = coords1(3,:) - minZ; %shifting z to be zero
        coords0 = coords0 + pos; % turn the x to be 130, shift the y
        coords1 = coords1 + pos; % turn the x to be 130, shift the y
        %gets the number of columns, i.e., the number of points around the boundary
        %This number is important for building the side walls of the extruded mesh
        n = size(coords0, 2);
        v0 = coords0'; %This represents the bottom face of the extruded letter.
        v1 = coords1'; %Same for the top face
        all_text_vertices{face_idx} = [v0; v1]; %storing both layers of vertices in your cell array
        faces_temp = []; %Initializes an empty list of triangle faces.
        %connecting side faces
        for m = 2:n-1
            faces_temp(end+1,:) = [vertex_offset+1, vertex_offset+m, vertex_offset+m+1];
            faces_temp(end+1,:) = [vertex_offset+n+1, vertex_offset+n+m, vertex_offset+n+m+1];
        end
        %extruding side walls
        for m = 1:n-1
            a = vertex_offset + m;
            b = vertex_offset + m + 1;
            c = vertex_offset + n + m + 1;
            d = vertex_offset + n + m;
            faces_temp(end+1,:) = [a b c];
            faces_temp(end+1,:) = [a c d];
        end
        %storing the triangulated faces and matrices
        all_text_faces{face_idx} = faces_temp;
        vertex_offset = vertex_offset + 2 * n;
        face_idx = face_idx + 1;
    end
%make the y ready for the next letter
    startY = startY + spacing;
end
%merging all the individual letters' geometry
all_text_vertices = vertcat(all_text_vertices{:});
all_text_faces = vertcat(all_text_faces{:});

%% Add Multiple Numbers
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

%% Add Custom Object on Top of Cube
% Assume user provides: f_solid, v_solid
% Position custom object centered on top of cube (e.g., 75, 100, 20)
custom_translation = [110, 20, 12]; 
rotation_xy = [cosd(90) -sind(90) 0; sind(90) cosd(90) 0; 0 0 1];
v_solid_rotated = (rotation_xy * V_solid')';
%moving 3D object
v_solid_translated = bsxfun(@plus, v_solid_rotated, custom_translation);
f_solid_shifted = F_solid + size(Solid_vertices, 1) + size( ...
    all_text_vertices, 1) + size(num_vertices, 1);

%% Combine All and Export
text_faces_shifted = all_text_faces + size(Solid_vertices, 1);
num_faces_shifted = num_faces + size(Solid_vertices, 1) + size(all_text_vertices, 1);
combined_V = [Solid_vertices; all_text_vertices; num_vertices; v_solid_translated];
combined_F = [Solid_faces; text_faces_shifted; num_faces_shifted; f_solid_shifted];
% 
TR = triangulation(combined_F, combined_V);
% stlwrite(TR, 'cube_with_text_number_and_custom_object.stl');
% disp('✅ Exported as cube_with_text_number_and_custom_object.stl');
% 
% %% Preview
% figure;
% patch('Faces', combined_faces, 'Vertices', combined_vertices, ...
%       'FaceColor', [0.8 0.3 1], 'EdgeColor', 'none');
% axis equal;
% camlight;
% lighting gouraud;
% title('Cube with 3D Text, Numbers, and Custom Object');
% xlabel('X'); ylabel('Y'); zlabel('Z');
end % the main function




 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 