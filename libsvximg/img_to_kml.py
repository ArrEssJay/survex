#!/usr/bin/python

from img import *
from UserDict import UserDict
import argparse, sys, math, cmath
from pyproj import Proj, Geod
from pykml.factory import KML_ElementMaker as KML
from lxml import etree
from vec3 import *
from math import *

#default config
epsg="epsg:28356" #GDA94 (UTM)
geod_sphere="GRS80"
#refstation = "ent.j13"
#ref_E = 0224154.45
#ref_N = 6255362.60
#ref_Z = 842.48
alt_offset = 300
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

def parse_3d(infile):
    "Parse a .3D file to lists of stations, cross sections and traverses"
    survey = ''
    pimg = img_open_survey(infile, survey)

    #print the .3d header
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

    #prime lists of stations, centreline traverses, tube traverses
    stations = {}
    traverses = []
    tubes = []

    #prime lists of stations in current tube/centreline traverses
    cXSTrav = []
    cLTrav = []


    #current traverse/tube state
    #cTrav is used to check that we didn't attempt to start a line before a move
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
                print("ERROR: not in a tube to end")
                #return 1
            tubes.append(cXSTrav)
            print("Added a tube with %d cross sections", len(cXSTrav))
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
                print("ERROR: Line before move")
                #return 1
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
    return traverses, tubes, stations


class XSect:
    "store cross section data"
    def __init__(self, l=None, r=None, u=None, d=None, station=None):
        self.l = l
        self.r = r
        self.u = u
        self.d = d
        self.station = station

class Tpos:
    INVALID = -1
    FIRST = 0
    INTER = 1
    LAST = 2

def ptocart(p, rs, rLon, rLat, geod):
    "Given a reference station with local co-ordinates rs at location rLon, rLat, calculate Lat and Lon for a station with local co-ordinates p"

    #Calculate x,y distance of point from reference
    dx = p.x - rs.x
    dy = p.y - rs.y

    #for altitude, add the delta from the reference station's elevation
    #to the reference station's altitude
    dz = p.z - rs.z
    sAlt = ref_Z + dz

    #convert cartesian dx, dy to polar r, phi
    phi = degrees(atan2(dx, dy))
    r = hypot(dx, dy)

    sLon, sLat, sBAz = geod.fwd(rLon, rLat, phi, r, radians=False)

    return sLon, sLat, (sAlt*alt_exag)+alt_offset

def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(description='Dump Contents of Survex .3d file as text')
    parser.add_argument('-i', '--input', help='Input 3D File', required=True)
    parser.add_argument('-r', '--ref_station', help='Station to georeference model to', required=True)
    parser.add_argument('-rs', '--ref_station', help='Station name to georeference model to', required=True)
    parser.add_argument('-re', '--ref_east', help='Reference station location Easting', required=True)
    parser.add_argument('-rn', '--ref_north', help='Reference station location Northing', required=True)
    parser.add_argument('-ra', '--ref_alt', help='Reference station altitude (m)', required=True)
    parser.add_argument('-d', '--datum', help='Datum string to supply to Proj for georeferncing. Defaults to "epsg:28356" (#GDA94 (UTM))', required=False)
    parser.add_argument('-g', '--datum', help='S', required=False)


    parser.add_argument('-o', '--output', help='Output File', required=False)
    parser.add_argument('-z', '--zip', help='zip contents to (.KMZ); automatically true for Collada output', required=False)
    parser.add_argument('-c', '--collada', help='output geometry in Collada .dae format', required=False)


    shell_args, unparsed_args = parser.parse_known_args()

    #print shell_args.input

    traverses, tubes, stations = parse_3d(shell_args.input)

    #locate the reference station
    refstation = shell_args.ref_station
    print "Reference Station Name: %s" %(refstation)
    if stations.has_key(refstation):
        print ("Reference station found")
        rs =stations.get(refstation)
        print "Reference Station @ x:%8.2f y:%8.2f z:%8.2f" % (rs.x, rs.y, rs.z)

    else:
        print("No reference station found. Can't continue")
        return 1

    #Calculate distances of all stations from reference station

    #calculate absolute position of reference station
    pProj = Proj(init=epsg)
    rLon, rLat = pProj(ref_E, ref_N,inverse=True)

    print "Reference Station Lat: %8.2f degrees, Lon: %8.2f degrees" % (rLat, rLon)

    #setup Geod
    geod = Geod(ellps=geod_sphere)

    #open KML
    #create KML object
    #kml_object = KML.name("Test")
    #fld = KML.Folder()
    # create a KML file skeleton
    surveystyle="survey_style"
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
                KML.Polystyle(
                    KML.color('bfff4e34'),
                    KML.colorMode('random'),
                ),
                id=surveystyle
            )


        )
    )



    cart_stations = {}
    kf = KML.Folder(KML.name("Stations"))

    for k,v in stations.iteritems():

        sLon,sLat,sAlt = ptocart(v, rs, rLon, rLat, geod)


        pm=KML.Placemark(
            KML.styleUrl('#{0}'.format(surveystyle)),

            KML.name(k),
            KML.description("Survex x,y,z = %8.2f %8.2f %8.2f" % (v.x, v.y, v.z)),
            KML.Point(
                       KML.coordinates("%8.8f,%8.8f,%8.8f" % (sLon,sLat,sAlt)),
                       KML.extrude(0),
                       KML.altitudeMode("absolute")
                    )
            )
        kf.append(pm)
    doc.Document.append(kf)


    #centreline
    kf = KML.Folder(KML.name("Centreline"))

    for tr, traverse in enumerate(traverses):
        print ("traverse %d of %d" % (tr, len(traverses)))
        pm =KML.Folder(
            KML.Name("Centreline"),
            KML.Placemark(
            KML.name("%s" % (tr)),
            KML.LineString(
                KML.coordinates(''.join(['%8.8f,%8.8f,%8.8f ' % (ptocart(leg, rs, rLon, rLat, geod)) for leg in traverse])
                                ),
                KML.altitudeMode("absolute")
                )
            ))
        kf.append(pm)

    doc.Document.append(kf)


    ######Now calculate the LRUD tubes
    up_v = vec3(0.0,0.0,1.0)
    last_right = vec3(0.0,0.0,0.0)
    U = [last_right,last_right,last_right,last_right]

    #setup KML folder
    kf = KML.Folder(KML.name("Tubes"))


    for tn,tube in enumerate(tubes):
        print "\n--------\nProcessing Tube: %d" %(tn)
        print "%d Segments" % len(tube)
        z_pitch_adjust = 0.0

        for vn, v_pres in enumerate(tube):

            right = vec3(0,0,0)
            up = vec3(0,0,0)
            tubepos = Tpos.INVALID

            #print right
            if vn == 0:
                #if this is the first segment we need to calculate the present vertex
                tubepos = Tpos.FIRST

                #why is XSect treated as a list here?
                #Calculate 'current' and 'next' vectors
                v_pres = v_pres[0]
                print "========\nVertex: %s" %(vn)
                print "Vertex pres: l: %8.2f r: %8.2f u: %8.2f d: %8.2f station: %s" % (v_pres.l, v_pres.r, v_pres.u, v_pres.d, v_pres.station)
                #get the X,Y,Z for this station
                v_pres_pos = stations[v_pres.station]
                v_pres_v = vec3(v_pres_pos.x, v_pres_pos.y, v_pres_pos.z)
                print "Vertex pres: %s" % (str(v_pres_v))

            #otherwise shift down and calculate the next segment
            else:
                #intermediary or last segment?
                if vn < len(tube)-1:
                    tubepos = Tpos.INTER
                else:
                    tubepos = Tpos.LAST

                v_prev_v = v_pres_v
                v_pres_v = v_next_v

                v_prev_pos = v_pres_pos
                v_pres_pos = v_next_pos

                #we still need the present lrud
                v_pres = v_pres[0]


            #if not the last segment, calculate next
            if tubepos == Tpos.FIRST or tubepos == Tpos.INTER:
                print("Need 'next' Vertex")
                v_next = tube[vn + 1]
                v_next = v_next[0]
                print "Vertex next: l: %8.2f r: %8.2f u: %8.2f d: %8.2f station: %s" % (v_next.l, v_next.r, v_next.u, v_next.d, v_next.station)
                v_next_pos = stations[v_next.station]
                v_next_v = vec3(v_next_pos.x, v_next_pos.y, v_next_pos.z)
                print "Vertex 'next' %s" % (str(v_next_v))

            #This vector rotation code is ported from the C++ implementation in gfxcore.cc
            #If this is the first segment we don't have the previous vector to project from
            if tubepos == Tpos.FIRST:
                print "First Segment"
                leg_v = vec3(v_next_v - v_pres_v)
                right = cross(leg_v,up_v)
                if abs(right) == 0:
                    right = last_right

                    up = up_v
                else:
                    last_right = right
                    up = up_v

            elif tubepos == Tpos.LAST:
                print "Last Segment"
                leg_v = vec3(v_pres_v - v_prev_v)
                right = cross(leg_v, up_v)


                if abs(right) == 0:
                    right = vec3(last_right.x, last_right.y, 0.0)

                    up = up_v
                else:
                    last_right = right
                    up = up_v


            elif tubepos == Tpos.INTER:
                print "Intermediary Segment"
                leg1_v = vec3(v_pres_v - v_prev_v)
                leg2_v = vec3(v_next_v - v_pres_v)

                r1 = cross(leg1_v, up_v)
                r2 = cross(leg2_v, up_v)


                r1 = r1.normalize()
                r2 = r2.normalize()


                right = vec3(r1 + r2)


                if abs(right) == 0:
                    #"mid-pitch case"
                    right = last_right

                if abs(r1) == 0:
                    n = leg1_v.normalize()
                    z_pitch_adjust = n.z
                    up = up_v
                    shift = 0
                    maxdotp = 0.0
                    right = right.normalize()
                    up = up.normalize()

                    vec = vec3(up - right)

                    for orient in range(0, 3):
                        tmp = vec3(U[orient] - v_prev_v)
                        tmp = tmp.normalize()
                        dotp = vec * tmp
                        if dotp > maxdotp:
                            maxdotp = dotp
                            shift=orient

                    if shift:
                        if shift !=2:

                            temp = U[0]
                            U[0] = U[shift]
                            U[shift] = U[2]
                            U[2] = U[shift^2]
                            U[shift^2] = temp
                        else:
                            U[0],U[2]=U[2],U[0]
                            U[1],U[3]=U[3],U[1]

                elif abs(r2) == 0:
                    n = leg2_v
                    n = n.normalize()
                    z_pitch_adjust = n.z
                    up = up_v

                else:
                    up = up_v

                last_right = right


            right = right.normalize()
            up = up.normalize()

            if z_pitch_adjust != 0:
                up = up + vec3(0, 0, abs(z_pitch_adjust))

            l = abs(v_pres.l)
            r = abs(v_pres.r)
            u = abs(v_pres.u)
            d = abs(v_pres.d)

            v = [(v_pres_v - (right * l) + (up * u)),
                    (v_pres_v + (right * r) + (up * u)),
                    (v_pres_v + (right * r) - (up * d)),
                    (v_pres_v - (right * l) - (up * d))
                    ]


            #now convert to lat,lon

            vll = [0,0,0,0]
            Ull = [0,0,0,0]
            for vit, vp in enumerate(v):
                tubept = img_point()
                tubept.x = vp.x
                tubept.y = vp.y
                tubept.z = vp.z

                sLon,sLat,sAlt = ptocart(tubept, rs, rLon, rLat, geod)
                vll[vit] =[sLon, sLat, sAlt]
            #print vll
            #print "vll: %s" % (str(vll))

            #print "U = %s" %(U)
            #print "v = %s" %(v)
            for Uit, Up in enumerate(U):
                tubept = img_point()
                tubept.x = Up.x
                tubept.y = Up.y
                tubept.z = Up.z

                sLon,sLat,sAlt = ptocart(tubept, rs, rLon, rLat, geod)
                Ull[Uit] =[sLon, sLat, sAlt]
            #print "Ull:  %s" % (str(Ull))
            if vn>0:
                for q in range(0,4):
                    #if q==0:
                    #    co = [vll[0], vll[1], Ull[1], Ull[0], vll[0]]
                    #elif q==1:
                    #    co = [vll[2], vll[3], Ull[3], Ull[2], vll[2]]
                    #elif q==2:
                    #    co = [vll[1], vll[2], Ull[2], Ull[1], vll[1]]
                    #elif q==3:
                    #    co = [vll[3], vll[0], Ull[0], Ull[3], vll[3]]

                    if q==0:
                        co = [vll[0], Ull[0], Ull[1], vll[1], vll[0]]
                    elif q==1:
                        co = [vll[2], Ull[2], Ull[3], vll[3], vll[2]]
                    elif q==2:
                        co = [vll[1], Ull[1], Ull[2], vll[2], vll[1]]
                    elif q==3:
                        co = [vll[3], Ull[3], Ull[0], vll[0], vll[3]]

                    pm = KML.Placemark(

                    KML.styleUrl('#{0}'.format(surveystyle)),
                    KML.name("tube: %s, vertex: %s, quad: %s" % (str(tn), str(vn), str(q))),
                    KML.Polygon(
                            KML.outerBoundaryIs(
                                KML.LinearRing(
                        KML.coordinates(''.join( ["%8.8f,%8.8f,%8.8f " % (leg[0], leg[1], leg[2]) for leg in co]))
                                    )),
                        KML.altitudeMode("absolute")
                        ))
                    kf.append(pm)

            if tubepos == Tpos.FIRST or tubepos == Tpos.LAST:
                #draw end cap
                cap = [vll[3], vll[2], vll[1], vll[0]]

                pm = KML.Placemark(

                    KML.styleUrl('#{0}'.format(surveystyle)),
                    KML.name("tube: %s end cap" % (str(tn))),
                    KML.Polygon(
                            KML.outerBoundaryIs(
                                KML.LinearRing(
                        KML.coordinates(''.join( ["%8.8f,%8.8f,%8.8f " % (leg[0], leg[1], leg[2]) for leg in cap]))
                                    )),
                        KML.altitudeMode("absolute")
                        ))
                kf.append(pm)


            doc.Document.append(kf)


#if (draw) {
##        if (segment > 0) {
 #       (this->*AddQuad)(v[0], v[1], U[1], U[0]);
 ##       (this->*AddQuad)(v[2], v[3], U[3], U[2]);
  #      (this->*AddQuad)(v[1], v[2], U[2], U[1]);
  #      (this->*AddQuad)(v[3], v[0], U[0], U[3]);
  #      }
#
 #       if (cover_end) {
  #      (this->*AddQuad)(v[3], v[2], v[1], v[0]);


     #assign after drawing



            U[0] = v[0]
            U[1] = v[1]
            U[2] = v[2]
            U[3] = v[3]




    outfile = file(__file__.rstrip('.py')+'.kml','w')
    outfile.write(etree.tostring(doc, pretty_print=True))



    #print tubes

    return 0

if __name__ == "__main__":
    sys.exit(main())



