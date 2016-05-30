readData;
%% Output setting for isabel_fld
max_x = 700;
max_y = 700;
unnormalized_field = isabel_fld;
field = normalizeField(isabel_fld);
field_name = 'isabel';
%% Output setting for isabel_fld
max_x = 300;
max_y = 300;
field_name = 'metsim1';
unnormalized_field = metsim1_fld;
field = normalizeField(metsim1_fld);

%% Plot vector fields
I_vectors = zeros([max_x max_y 3]);
I_vectors = plotVectorField(field, I_vectors, 10, 10, 6);
%% Plot single stream line with different step_size
step_size = 0.5;
start_x = 500; start_y = 500;
max_length = 2000;
COLOR_1 = [0 1 1 1]; COLOR_2 = [80/255 1 1 1];  COLOR_3 = [120/255 1 1 1]; 
fieldline = forwardEuler( start_x, start_y, max_x, max_y, step_size, max_length, field);
I1 = bitmapplot(fieldline.x, fieldline.y, I_vectors, struct('LineWidth',1,'Color',COLOR_1));
figure, imshow(hsv2rgb(I1));
fieldline = rungeKutta( start_x, start_y, max_x, max_y, step_size, max_length, field);
I2 = bitmapplot(fieldline.x, fieldline.y, I_vectors, struct('LineWidth',1,'Color',COLOR_2));
figure, imshow(hsv2rgb(I2));
fieldline = adaptiveRungeKutta( start_x, start_y, max_x, max_y, step_size, max_length, field, false, 0.1);
I3 = bitmapplot(fieldline.x, fieldline.y, I_vectors, struct('LineWidth',1,'Color',COLOR_3));
figure, imshow(hsv2rgb(I3));

step_size = 10;
COLOR_1 = [0 1 1 1]; COLOR_2 = [80/255 1 1 1];  COLOR_3 = [120/255 1 1 1]; 
fieldline = forwardEuler( start_x, start_y, max_x, max_y, step_size, max_length, field);
I1 = bitmapplot(fieldline.x, fieldline.y, I_vectors, struct('LineWidth',1,'Color',COLOR_1));
figure, imshow(hsv2rgb(I1));
fieldline = rungeKutta( start_x, start_y, max_x, max_y, step_size, max_length, field);
I2 = bitmapplot(fieldline.x, fieldline.y, I_vectors, struct('LineWidth',1,'Color',COLOR_2));
figure, imshow(hsv2rgb(I2));
fieldline = adaptiveRungeKutta( start_x, start_y, max_x, max_y, step_size, max_length, field, false, 0.1);
I3 = bitmapplot(fieldline.x, fieldline.y, I_vectors, struct('LineWidth',1,'Color',COLOR_3));
figure, imshow(hsv2rgb(I3));
%% Plot stream lines
step_size = 0.5;
%Grid
I_streamLineEuler = plotStreamLinesGrid(max_x, max_y, step_size, field, false, @forwardEuler);
%I_streamLineRK = plotStreamLinesGrid(max_x, max_y, step_size, field, false, @rungeKutta);
%I_streamLineAdaptiveRK = plotStreamLinesGrid(max_x, max_y, step_size, field, false, @adaptiveRungeKutta);
%Poisson Disc
I_streamLinePoission = plotStreamLinesPoisson(max_x, max_y, step_size, field, false, @adaptiveRungeKutta);
%Magnitude of velocity as degree of interest
p = zeros([max_x max_y 1]);
p = calProperty(p, unnormalized_field, 'magnitude_Velocity');
plotStreamLinesProp(p, max_x, max_y, step_size, field, false, @adaptiveRungeKutta);
%Enstrophy as degree of interest
p = zeros([max_x max_y 1]);
p = calProperty(p, unnormalized_field, 'enstrophy');
plotStreamLinesProp(p, max_x, max_y, step_size, field, false, @adaptiveRungeKutta);
%Topology-based
%Since the topoDegree cannot always correctly identify critical points,
%here I manually write the correct critical point with type to demonstrate
%the topo-based method
TopoDegree = zeros([max_x max_y 1]);
if strcmp(field_name, 'isabel')
    TopoDegree(90, 437) = 1; %saddle
    TopoDegree(538, 397) = 2; %spiral
    I_streamLineTopo = plotStreamLinesTopobased(TopoDegree, max_x, max_y, step_size, field, false, @adaptiveRungeKutta);
else
    TopoDegree(150, 150) = 2; %spiral
    I_streamLineTopo = plotStreamLinesTopobased(TopoDegree, max_x, max_y, step_size, field, false, @adaptiveRungeKutta);

end
%Maybe Better when merge with Poission
I = max(I_streamLineTopo, I_streamLinePoission);
figure, imshow(hsv2rgb(I));
%The evenly spaced streamline method needs careful implementation to be
%feasible, I ignore it here.
%% LIC
%Step_size is fixed at 0.5
input_texture = randi([0, 255], max_x, max_y)/255;
output = lineInterConv(input_texture, field, @forwardEuler);
figure, imshow(output.I);
output = lineInterConv(input_texture, field, @rungeKutta); %Really slow
figure, imshow(output.I);
output = lineInterConv(input_texture, field, @adaptiveRungeKutta);
figure, imshow(output.I);
%% Illuminated Stream Lines;
I_illuminated = plotStreamLinesGrid(max_x, max_y, step_size, field, true, @adaptiveRungeKutta);

%% Animation Streamlines;
animateStreamLine(max_x, max_y, field);
%% Topological Degree
TopoDegree = zeros([max_x max_y 1]);
TopoDegree = calProperty(TopoDegree, field, 'topoDegree');
figure, imshow(TopoDegree);
%% Classify Critical points using Trace-Determinant Plane
TypeCriticalPoint = zeros([max_x max_y 3]);
color_map = [[0 0 0]; [0 0 1]; [0 1 0]; [0 1 1]; [1 0 0]; [1 0 1]];
idx = find(TopoDegree>0);

for i = 1:length(idx)
    [x, y] = ind2sub(size(TopoDegree), idx(i));
    type = classifyCritical(x, y, max_x, max_y, unnormalized_field);
    TypeCriticalPoint(x, y, :) = color_map(type+1, :);
    TypeCriticalPoint = bitmaptext(num2str(type), TypeCriticalPoint, [x, y]);
end
figure, imshow(TypeCriticalPoint);

%% Enstrophy
Enstrophy = zeros([max_x max_y 1]);
Enstrophy = calProperty(Enstrophy, unnormalized_field, 'enstrophy');
figure, imshow(sqrt(Enstrophy));
%% Q criteria
Q_criteria = zeros([max_x max_y 1]);
Q_criteria = calProperty(Q_criteria, unnormalized_field, 'Q_criteria');
figure, imshow(Q_criteria);
