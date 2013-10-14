#!/usr/bin/python

from img import *
import argparse, sys

def main(argv=None):
    global dryrun, thexec, adbexec, waitfordevice
    if argv is None:
        argv = sys.argv

    print argv

    parser = argparse.ArgumentParser(description='Dump Contents of Survex .3d file as text')
    parser.add_argument('-i', '--input', help='Input 3D File', required=True)

    shell_args, unparsed_args = parser.parse_known_args()

    print shell_args.input
    survey = ''
    pimg = img_open_survey(shell_args.input, survey)

    if(pimg):
        print('File: '+ shell_args.input)
        print('------\nHeader\n------')
        print('Title: ' + str(pimg.title))
        print('.3d Version: ' + str(pimg.version))
        print('Datestamp: ' + str(pimg.datestamp))
        print('Level Separator: ' + str(pimg.separator))
        if (pimg.flags != 0):
            print('Flags: Extended Elevation')
    else:
        print('Failed to open: ' + shell_args.input)
        return 1

    stations = {}

    item = 0
    count_move = count_line = count_label = count_xsect = count_einfo = 0
    p_xyz = p_lrud = p_label = p_flags = p_einfo = p_style = 0
    x_max = y_max = z_max = 0
    x_min = y_min = z_min = 0
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
            count_einfo += 1
        elif (item == img_XSECT):
            print("XSECT")
            p_lrud = p_label = p_flags = p_style = 1
            count_xsect += 1
        elif (item == img_XSECT_END):
            print("XSECT_END")
            p_label = p_flags = 1
        elif (item == img_MOVE):
            print("MOVE")
            p_xyz = p_label = p_flags = 1
            count_move += 1
        elif (item == img_LINE):
            print("LINE")
            p_xyz = p_label = p_flags = p_style = p_lrud = 1
            count_line += 1
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

            #update x/y/z_max
            #Set an initial value
            if(count_label == 0):
                x_max = x_min = pt.x
                y_max = y_min = pt.y
                z_max = z_min = pt.z
            else:
                if (pt.x > x_max):
                    x_max = pt.x
                if (pt.y > y_max):
                    y_max = pt.y
                if (pt.z > z_max):
                    z_max = pt.z

                if (pt.x < x_min):
                    x_min = pt.x
                if (pt.y < y_min):
                    y_min = pt.y
                if (pt.z < z_min):
                    z_min = pt.z
            count_label += 1

        if (p_xyz):
            print("     X: "+str(p.x) +" Y: "+str(p.y) +" Z: "+str(p.z))
            p_xyz = 0

        if (p_lrud):
            print("     LRUD: l:"+str(pimg.l)+" r:"+str(pimg.r)+" u:"+str(pimg.u)+" d:"+str(pimg.d))
            p_lrud = 0

        if (p_einfo):
            print("     ERROR: # Legs: "+str(pimg.n_legs)+" Length: "+str(pimg.length))
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


    


    print("Station Labels->XYZ:")
    for k, pt in stations.iteritems():
        print(k+" -> X: "+str(pt.x)+" Y: "+str(pt.y)+" Z: "+str(pt.z))

    print ("------\nStats\n-----")
    print("Total LABEL: "+str(count_label))
    print("Total LINE: "+str(count_line))
    print("Total XSECT: "+str(count_xsect))
    print("Total MOVE: "+str(count_move))
    print("Total ERROR_INFO: "+str(count_einfo))
    print('-----')
    #Calculate X/Y/Z ranges
    print ("Max: x: "+str(x_max)+" y: "+str(y_max)+" z: "+str(z_max))
    print ("Min: x: "+str(x_min)+" y: "+str(y_min)+" z: "+str(z_min))
    print ("Range: x: "+str(x_max-x_min)+" y: "+str(y_max-y_min)+" z: "+str(z_max-z_min))

    return 0

if __name__ == "__main__":
    sys.exit(main())
