function [center, radius] = fitcenter(inputArg1)
    v = VideoReader(inputArg1);
    frame = read(v, 1);
    gray_frame = rgb2gray(frame);
    
    imshow(edge(gray_frame, 'sobel', 0.0125))
    c = drawcircle;
    pause;
    center = c.Center;
    radius = c.Radius;
end