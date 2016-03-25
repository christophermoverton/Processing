## to be used in blender
import bpy
import csv
import os
import sys

obj_list = ["CurveObj", "CurveObj2", "CurveObj3", "CurveObj4","CurveObj5"]
ob_curve = bpy.data.objects["CurveObj3"]
cu = ob_curve.data
points = []
lpoints = []
rpoints = []
for obj in obj_list:
    ob_curve = bpy.data.objects[obj]
    cu = ob_curve.data
    for spline in cu.splines:
        for point in spline.bezier_points:
            points.append(tuple(point.co))
            lpoints.append(tuple(point.handle_left))
            rpoints.append(tuple(point.handle_right))
print(points)
print(lpoints)
print(rpoints)
print(bpy.context.space_data.text.filepath)
filepath = bpy.context.space_data.text.filepath.split('beziercurvedataread2.py')[0]
print(filepath)
with open(filepath+'bezierpoints.csv', 'w') as csvfile:
    fieldnames = ['x', 'y', 'z', 'cplx', 'cply', 'cplz',
                  'cprx', 'cpry', 'cprz']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for i,point in enumerate(points):
        writer.writerow({'x': point[0], 'y': point[1], 'z': point[2],
                         'cplx': lpoints[i][0], 'cply': lpoints[i][1], 
                         'cplz':lpoints[i][2], 'cprx': rpoints[i][0],
                         'cpry':rpoints[i][1], 'cprz': rpoints[i][2]})