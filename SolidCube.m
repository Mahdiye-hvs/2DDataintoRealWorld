% Define the 8 vertices of the cube
vertices = [0 0 0;
            100 0 0;
            100 150 0;
            0 150 0;
            0 0 5;
            100 0 5;
            100 150 5;
            0 150 5];

faces = delaunay(vertices(:, 1), vertices(:, 2)); %faces trangles

[Solid_faces, Solid_vertices] = surf2solid(faces, vertices, 'elevation', 10); %it should have the triangles faces to solidify

text(topCenter(1), topCenter(2), topCenter(3), 'Top Face', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'Color', 'k')

% Plot the cube
figure;
patch('Faces', Solid_faces, 'Vertices', Solid_vertices,'FaceColor', 'cyan', 'EdgeColor', 'none');
axis vis3d;
view(3);
camlight;
lighting gouraud;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Solid Volume');

% Create a triangulation object from decimated solid 
tri_object = triangulation(Solid_faces, Solid_vertices); %needed for exporting data to stl format

% Write the decimated triangulation object to an STL file
stlwrite(tri_object, 'LC×LC_StlFormat.stl');



