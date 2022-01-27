v = linspace(-20, 20, 101);
figure()
hold on
plot(WALL.RWALL_list, WALL.ZWALL_list, '-k'); % 容器壁 VacuumVesselMeshPoints
contour(REF.R, REF.Z, REF.Flux*1000,v, 'LineColor', 'm');