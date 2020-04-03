function SLarray = calcSL(marker1,marker2,indHS1)

% Find SL as distance between two Heel markers at one HS event
% Look at AP direction of markers only

SLarray = marker1(indHS1) - marker2(indHS1);