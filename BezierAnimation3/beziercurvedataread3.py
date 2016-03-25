## to be used in blender
import bpy
import csv
import os
import sys
Scene_Name = bpy.context.scene.name
##obj_list = ["CurveObj", "CurveObj2", "CurveObj3", "CurveObj4","CurveObj5"]
obj_list = ["CurveObj"]
hook_list = ["hook1", "hook2", "hook3", "hook4", "hook5"]
##ob_curve = bpy.data.objects["CurveObj3"]
##cu = ob_curve.data
points_animation = []
points = {}
lpoints = {}
rpoints = {}
for frame in range(1, 120):
    bpy.data.scenes[Scene_Name].frame_set(frame)
    for i, obj in enumerate(obj_list):
        ob_curve = bpy.data.objects[obj]
        cu = ob_curve.data
        for spline in cu.splines:
            maxPts = len(spline.bezier_points)
            for j,point in enumerate(spline.bezier_points):
                ob_hook = bpy.data.objects[hook_list[j]]
                points[len(points)] = {'pos': tuple(ob_hook.location.co), 
                                        'curveID': i, 'frameID': frame,
                                        'pointID': j, 'maxPts': maxPts}
                lpoints[len(lpoints)] = {'pos':tuple(point.handle_left),
                                        'curveID': i, 'frameID': frame,
                                        'pointID': j}
                rpoints[len(rpoints)] = {'pos':tuple(point.handle_right),
                                        'curveID': i, 'frameID': frame,
                                         'pointID': j}
print(points)
print(lpoints)
print(rpoints)
print(bpy.context.space_data.text.filepath)
filepath = bpy.context.space_data.text.filepath.split('beziercurvedataread2.py')[0]
print(filepath)
filepath = '/home/strangequark/CProcessing/BezierAnimation3/'
with open(filepath+'bezierpoints.csv', 'w') as csvfile:
    fieldnames = ['x', 'y', 'z', 'cplx', 'cply', 'cplz',
                  'cprx', 'cpry', 'cprz', 'pointID','frameID', 'curveID',
                  'maxPts']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for point in points:
        writer.writerow({'x': points[point]['pos'][0], 
                         'y': points[point]['pos'][1], 
                         'z': points[point]['pos'][2],
                         'cplx': lpoints[point]['pos'][0], 
                         'cply': lpoints[point]['pos'][1], 
                         'cplz':lpoints[point]['pos'][2], 
                         'cprx': rpoints[point]['pos'][0],
                         'cpry':rpoints[point]['pos'][1], 
                         'cprz': rpoints[point]['pos'][2],
                         'pointID': points[point]['pointID'],
                         'frameID': points[point]['frameID'],
                         'curveID': points[point]['curveID'],
                         'maxPts': points[point]['maxPts']})