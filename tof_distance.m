function X = tof_distance(x)
% Get the distance between sample and mass spectrometer from the distance
% "measured on the outside".

p1 = 0.91319; p2 = -33.8929; 
X = p1*x + p2;