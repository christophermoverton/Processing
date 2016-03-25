import bpy
obj_list = ["CurveObj"]
#obj = bpy.data.objects.new("link", curve)

#m0 = obj.modifiers.new("alpha", 'HOOK')
for i, obj in enumerate(obj_list):
    ob_curve = bpy.data.objects[obj]
    cu = ob_curve.data
    for spline in cu.splines:
        modnames = ['hook.'+str(k).zfill(3) for k in range(len(spline.bezier_points))]
        lmodnames = ['lhook.'+str(k).zfill(3) for k in range(len(spline.bezier_points))]
        rmodnames = ['rhook.'+str(k).zfill(3) for k in range(len(spline.bezier_points))]    

#        spline.bezier_points[3].select_control_point = True
#        spline.bezier_points[3].select_left_handle = True
        print(modnames)
#    for i, modname in enumerate(modnames):
#        bpy.ops.object.add(type='EMPTY')
#        obj = bpy.data.objects["Empty"]
#        print(obj)
#        obj.name = modname
#  
#        obj.modifiers.new(modname, type='HOOK')
#        bpy.ops.object.add(type='EMPTY')
#        obj = bpy.data.objects["Empty"]
#        print(obj)
#        obj.name = lmodnames[i]
#  
#        obj.modifiers.new(rmodnames[i], type='HOOK') 
#        bpy.ops.object.add(type='EMPTY')
#        obj = bpy.data.objects["Empty"]
#        print(obj)
#        obj.name = rmodnames[i]
#  
#        obj.modifiers.new(rmodnames[i], type='HOOK') 
    ##bpy.ops.object.mode_set(mode = 'EDIT')
#    obj_curve.active = True
    bpy.ops.object.select_all(action='DESELECT') 
    ob_curve = bpy.data.objects["CurveObj"]
    ob_curve.select = True
    bpy.context.scene.objects.active = ob_curve
    bpy.ops.object.mode_set(mode='EDIT')
    for spline in cu.splines:
        for i, point in enumerate(spline.bezier_points):
            
            point.select_control_point = True
#            obj = bpy.data.objects[modnames[i]]
#            obj.location = point.co
            ##bpy.ops.object.hook_assign(modifier=modnames[i])
            bpy.ops.object.hook_add_newob()
            obj = bpy.data.objects["Empty"]
            print(obj)
            obj.name = modnames[i]
            point.select_control_point = False
            point.select_left_handle = True
#            obj = bpy.data.objects[lmodnames[i]]
#            obj.location = point.handle_left
#            bpy.ops.object.hook_assign(modifier=lmodnames[i])
            bpy.ops.object.hook_add_newob()
            obj = bpy.data.objects["Empty"]
            print(obj)
            obj.name = lmodnames[i]
            point.select_left_handle = False
            point.select_right_handle = True
#            obj = bpy.data.objects[rmodnames[i]]
#            obj.location = point.handle_right
#            bpy.ops.object.hook_assign(modifier=rmodnames[i])
            bpy.ops.object.hook_add_newob()
            obj = bpy.data.objects["Empty"]
            print(obj)
            obj.name = rmodnames[i]
            point.select_right_handle = False