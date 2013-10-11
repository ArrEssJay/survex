#!/usr/bin/python

from img import *

fnm_3d = '/Users/rob/Dropbox/Personal/caving/cavesstuff/phil_3d_sample/mammoth.3d'

survey = ''
pimg = img_open_survey(fnm_3d, survey)

if(pimg):
    print('File: '+ fnm_3d)
    print('------\nHeader\n------')
    print('Title: ' + str(pimg.title))
    print('.3d Version: ' + str(pimg.version))
    print('Datestamp: ' + str(pimg.datestamp))
    print('Level Separator: ' + str(pimg.separator))
    if (pimg.flags != 0):
        print('Flags: Extended Elevation')
else:
	print('Failed to open: '+fnm_3d)

stations = {}

item = 0
count = 0
p_xyz = p_lrud = p_label = p_flags = p_einfo = p_style = 0
p = img_point()

while(item != img_STOP):
    item = img_read_item(pimg, p);

    if (item == img_BAD):
        print("BAD")
    elif (item == img_STOP):
        print("STOP")
    elif (item == img_ERROR_INFO):
        print("ERROR_INFO")
        p_einfo = p_pxy = p_style = 1
    elif (item == img_XSECT):
        print("XSECT")
        p_lrud = p_label = p_flags = p_style = 1   
    elif (item == img_XSECT_END):
        print("XSECT_END")
        p_label = p_flags = 1 
    elif (item == img_MOVE):
        print("MOVE")
        p_xyz = p_label = p_flags = 1
    elif (item == img_LINE):
        print("LINE")
        p_xyz = p_label = p_flags = p_style = 1
    elif (item == img_LABEL):
        print("LABEL")
    
        #Can't just assign pt to p as it will change on next read_item()
        #Nor will deepcopy work without changes to the interface
        #Must be a more elgant way to do this
        pt = img_point()
        pt.x = p.x
        pt.y = p.y
        pt.z = p.z

        stations.update({pimg.label:pt})
        p_xyz = p_label = p_flags = p_style = 1

    if (p_xyz):
        print("     X: "+str(p.x) +" Y: "+str(p.y) +" Z: "+str(p.z))
        p_xyz = 0
        
    if (p_lrud):
        print("     LRUD: l:"+str(pimg.l)+" r:"+str(pimg.r)+" u:"+str(pimg.u)+" d:"+str(pimg.d))
        p_lrud = 0

    if (p_einfo):
        print("     ERROR: E: "+str(pimg.E)+" H: "+str(pimg.H)+" V: "+str(pimg.V))
        p_einfo = 0
    
    if (p_label):
        labels = pimg.label.split(pimg.separator)
        if (len(labels) >= 1 ):
            print("     LEVEL: " + labels[0])
        if (len(labels) == 2 ):
            print("     STATION: " + labels[1])
        p_label = 0
    
    if (p_flags and (pimg.flags != 0)):
        print("     FLAGS("+hex(pimg.flags)+"):"),
        if (item == img_LINE):
            if(pimg.flags & img_FLAG_SURFACE):
                print("LINE_FLAG_SURFACE"),
            if(pimg.flags & img_FLAG_DUPLICATE):
                print("LINE_FLAG_DUPLICATE"),
            if(pimg.flags & img_FLAG_DUPLICATE):
                print("LINE_FLAG_SPLAY"),
        elif (item == img_LABEL):
            if(pimg.flags & img_SFLAG_SURFACE):
                print("STATION_FLAG_SURFACE"),
            if(pimg.flags & img_SFLAG_UNDERGROUND):
                print("STATION_FLAG_UNDERGROUND"),
            if(pimg.flags & img_SFLAG_ENTRANCE):
                print("STATION_FLAG_ENTRANCE"),
            if(pimg.flags & img_SFLAG_EXPORTED):
                print("STATION_FLAG_EXPORTED"),
            if(pimg.flags & img_SFLAG_FIXED):
                print("STATION_FLAG_FIXED"),
            if(pimg.flags & img_SFLAG_ANON):
                print("STATION_FLAG_ANON"),
            if(pimg.flags & img_SFLAG_WALL):
                print("STATION_FLAG_WALL"),
        elif ((item == img_XSECT) or (item == img_XSECT_END)):
            if(pimg.flags & img_XFLAG_END):
                print("XSECT_FLAG_END"),
        else:
            print("UNKNOWN FLAG: "+hex(pimg.flags)),
        print('')
        p_flags = 0
    
    if(p_style and (pimg.style != img_STYLE_UNKNOWN)):
        print("     STYLE:"),
        if(pimg.style == img_STYLE_NORMAL):
            print("NORMAL"),
        elif(pimg.style == img_STYLE_DIVING):
            print("DIVING"),
        elif(pimg.style == img_STYLE_CARTESIAN):
            print("CARTESIAN"),
        elif(pimg.style == img_STYLE_CYLPOLAR):
            print("CYLPOLAR"),
        elif(pimg.style == img_STYLE_NOSURVEY):
            print("NOSURVEY"),
        print('')
        p_style = 0

    count += 1

print("Total Points: "+str(count))

print("Station Labels->XYZ:")
for k, pt in stations.iteritems():
    print(k+" -> X: "+str(pt.x)+" Y: "+str(pt.y)+" Z: "+str(pt.z))


