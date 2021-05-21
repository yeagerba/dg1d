% Define a Mesh class
classdef MeshClass
  properties
    Points
    Connectivity
    dx
    dx_min
    Nelems
  end
  methods
    function MC = MeshClass(points, connectivity)
      MC.Points = points;
      MC.Connectivity = connectivity;
      MC.dx = points(2:end) - points(1:end-1);
      MC.dx_min = min(MC.dx);
      MC.Nelems = length(points) - 1;
    end
  end
end
