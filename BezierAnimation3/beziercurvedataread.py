## to be used in blender
import bpy
ob_curve = bpy.data.objects["BezierCurve"]
cu = ob_curve.data
points = []
lpoints = []
rpoints = []
for spline in cu.splines:
    for point in spline.bezier_points:
        points.append(tuple(point.co))
        lpoints.append(tuple(point.handle_left))
        rpoints.append(tuple(point.handle_right))
print(points)
print(lpoints)
print(rpoints)