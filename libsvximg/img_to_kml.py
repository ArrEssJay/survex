#!/usr/bin/python

from img import *
from UserDict import UserDict
import argparse, sys, math, cmath
from pyproj import Proj, Geod
from pykml.factory import KML_ElementMaker as KML
from lxml import etree


#config
epsg="epsg:28356" #GDA94 (UTM)
refstation = "ent.j13"
ref_E = 0224154.45
ref_N = 6255362.60
ref_Z = 842.48
alt_offset = 200
alt_exag = 1

def copy_point(p):
    #Can't just assign pt to p as it will change on next read_item()
    #Nor will deepcopy work without changes to the interface
    #Must be a more elgant way to do this
    pt = img_point()
    pt.x = p.x
    pt.y = p.y
    pt.z = p.z
    return pt

class XSect(UserDict):
    "store cross section data"
    def __init__(self, l=None, r=None, u=None, d=None, station=None):
        UserDict.__init__(self)
        self["l"] = l
        self["r"] = r
        self["u"] = u
        self["d"] = d
        self["station"] = station

def ptocart(p, rs, rLon, rLat, geod):
    dx = p.x - rs.x
    dy = p.y - rs.y

    #for altitude, calculate the delta from the reference station
    #then add to the value at the reference mark
    dz = p.z - rs.z
    sAlt = ref_Z + dz
    print 'dx=%8.2f, dy=%8.2f, dz=%8.2f' % (dx,dy,dz)

    #avoid divide by zero
    if dx == 0:
        phi = 0
    else:
        #wrt x plane
        phi = 90 + math.degrees(math.atan(dy/dx))


    r = math.sqrt(pow(dx,2) + pow(dy,2))

    print"r=%8.2f, phi=%8.2f" % (r,phi)

    sLon, sLat, sBAz = geod.fwd(rLon, rLat, phi, r, radians=False)
    print"sLon = %8.2f, sLat = %8.2f, sBAz = %8.2f, sAlt = %8.2f" % (sLon, sLat, sBAz, sAlt)

    return sLon, sLat, (sAlt*alt_exag)+alt_offset

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

    #Parse the file and get lists of traverses


    stations = {}
    traverses = []
    tubes = []

    cXSTrav = []
    cLTrav = []


    #current traverse/tube state
    #cTrav is only to check if we did a line before a move
    cTrav = cTube = False

    p = img_point()

    item = 0
    while(item != img_STOP):
        item = img_read_item(pimg, p);

        if (item == img_BAD):
            print("BAD")
        elif (item == img_STOP):
            print("STOP")
        #elif (item == img_ERROR_INFO):

        elif (item == img_XSECT):
            xS = XSect(pimg.l, pimg.r, pimg.u, pimg.d, pimg.label)

            #clear the previous tube
            if(cTube == False):
                cXSTrav = []
            cXSTrav.append([xS])
            cTube = True

        elif (item == img_XSECT_END):
            if(cTube == False):
                print("Error: not in a tube to end")
                return 1
            tubes.append(cXSTrav)
            #print("Added a tube with %d cross sections", len(cXSTrav))
            cTube = False

        elif (item == img_MOVE):
            if(cTrav):
                #we already drew a line so this is the start of a second
                #move the previous traverse to the traverses list
                traverses.append(cLTrav)

            #start a new traverse
            cLTrav = []
            pt = copy_point(p)
            cLTrav.append(pt)
            cTrav = True

        elif (item == img_LINE):
            if(cTrav == False):
                print("Line before move")
                return 1
            pt = copy_point(p)
            cLTrav.append(pt)

        elif (item == img_LABEL):
            pt = copy_point(p)
            stations.update({pimg.label:pt})

    print "Total: %d Traverses" % (len(traverses))
    #for i1, traverse in enumerate(traverses):
        #print "    #%d : %d legs" % (i1, len(traverse))
        #for i2, leg in enumerate(traverse):
            #print "Leg: %d, : x:%8.2f, y:%8.2f, z:%8.2f" % (i2, leg.x, leg.y, leg.z)

    print "Total: %d Tubes" % (len(tubes))
    #for i1, tube in enumerate(tubes):
        #print "    #%d : %d legs" % (i1, len(tube))
        #for i2, pxs in enumerate(tube):
            #print "Leg: %d, : l:%8.2f, r:%8.2f, u:%8.2f, d:%8.2f" % (i2, pxs.l, pxs.r, pxs.u, pxs.d)
            #print pxs

    print "Total: %d Stations" %(len(stations))
    #print stations.items()
    #for k,v in stations.items():
    #    print ("%s : %8.2f %8.2f %8.2f") % (k,v.x, v.y, v.z)



        #find the reference station

    print "Refstation: %s" %(refstation)
    if stations.has_key(refstation):
        print ("Reference station found")
        rs =stations.get(refstation)
        print "Reference Station @ x:%8.2f y:%8.2f z:%8.2f" % (rs.x, rs.y, rs.z)

    else:
        print("No reference station found. Can't continue")
        return 1

    #Calculate distances of all stations from reference station

    #calculate absolute position of reference station
    p = Proj(init=epsg)
    rLon, rLat = p(ref_E, ref_N,inverse=True)

    print "Reference Station Lat: %8.2f degrees, Lon: %8.2f degrees" % (rLat, rLon)

    #setup Geod
    geod = Geod(ellps='GRS80')

    #open KML
    #create KML object
    #kml_object = KML.name("Test")
    #fld = KML.Folder()
    # create a KML file skeleton
    stylename="survey_style"
    doc = KML.kml(
        KML.Document(
            KML.Name(str(pimg.title)),
            KML.Style(
                KML.IconStyle(
                    KML.scale(0.4),
                    KML.Icon(
                        KML.href("http://maps.google.com/mapfiles/kml/pal4/icon49.png")
                    ),
                ),

                KML.LabelStyle(
                    KML.scale(0.6)
                ),
                id=stylename
            )


        )
    )

    cart_stations = {}
    for k,v in stations.iteritems():

        sLon,sLat,sAlt = ptocart(v, rs, rLon, rLat, geod)

        pm = KML.Placemark(
            KML.styleUrl('#{0}'.format(stylename)),

            KML.name(k),
            KML.description("Survex x,y,z = %8.2f %8.2f %8.2f" % (v.x, v.y, v.z)),
            KML.Point(
                       KML.coordinates("%8.8f,%8.8f,%8.8f" % (sLon,sLat,sAlt)),
                       KML.extrude(0),
                       KML.altitudeMode("absolute")
                    )
            )
        doc.Document.append(pm)


    #centreline
    for tr, traverse in enumerate(traverses):
        print ("traverse %d of %d" % (tr, len(traverses)))
        pm = KML.Placemark(
            KML.name("%s" % (tr)),
            KML.LineString(
                KML.coordinates(''.join(['%8.8f,%8.8f,%8.8f ' % (ptocart(leg, rs, rLon, rLat, geod)) for leg in traverse])
                                ),
                KML.altitudeMode("absolute")
                )
            )
        doc.Document.append(pm)


    outfile = file(__file__.rstrip('.py')+'.kml','w')
    outfile.write(etree.tostring(doc, pretty_print=True))






    #print tubes

    return 0

if __name__ == "__main__":
    sys.exit(main())



