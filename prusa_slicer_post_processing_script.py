"""
This script generates Overhangs by stringing together Arcs, allowing successful fdm-3d-printing of large 90 deg overhangs!
The genius Idea is emerged from Steven McCulloch, who coded a demonstration and the basic mechanics: https://github.com/stmcculloch/arc-overhang
This python script builds up on that and offers a convinient way to integrate the ArcOverhangs into an existing gcode-file.

HOW TO USE: 
Option A) open your system console and type 'python ' followed by the path to this script and the path of the gcode file. Will overwrite the file.
Option B) open PrusaSlicer, go to print-settings-tab->output-options. Locate the window for post-processing-script. 
    In that window enter: full path to your python exe,emtyspace, full path to this script.
    If the python path contains any empty spaces, mask them as described here: https://manual.slic3r.org/advanced/post-processing
=>PrusaSlicer will execute the script after the export of the Gcode, therefore the view in the window wont change. Open the finished gcode file to see the results.

If you want to change generation settings: Scroll to 'Parameter' section. Settings from PrusaSlicer will be extracted automaticly from the gcode.

Requirements:
Python 3.5+ and the librarys: shapely 1.8+, numpy 1.2+, matplotlib for debugging

Limitations:
-The Arc-Generation is not hole-tolerant yet.
-All Usecases with multiple spots for arc generation are not testet yet and might contain bugs.
-The Arcs are extruded very thick, so the layer will be 0.1-0.5mm thicker (dependend on nozzle dia) than expected=>if precision needed make test prints to counter this effect.

Notes:
This code is a little messy. Usually I would devide it into multiple files, but that would compromise the ease of use.
Therefore I divided the code into sections, marked with ###
Feel free to give it some refactoring and add more functionalities!

Future possible extensions:
-slow down all commands 5mm above arcs, to prevent warping, slow down the bridging-speed of the perimeter, if necessary
-hole tolerant
-print the small arcs even slower with dynamic speed adjustment

Known issues:
-pointsPerCircle>80 give weird results
-maxDistanceFromPerimeter >=2*perimeterwidth=weird result.
-avoid using the code multiple times onto the same gcode, as errors might accumulate.
"""
#!/usr/bin/python
import sys
import os
from shapely import Point, Polygon, LineString, GeometryCollection, MultiLineString, MultiPolygon
from shapely.ops import nearest_points
from shapely.ops import linemerge
import matplotlib.pyplot as plt
import numpy as np
from ast import literal_eval
import warnings

########## Parameters  - adjust values here as needed ##########
def makeFullSettingDict(gCodeSettingDict:dict) -> dict: 
    """Merge Two Dictionarys and set some keys/values explicitly"""
    #add keys below, when you want to add a setting or overwrite an existing key from the config in the gcode. For some parameters the default settings form the gcode can be copied via gCodeSettingDict.get(key)
    AddManualSettingsDict={
        "notPerformPerimeterCheck":False, #Usual gcode:external perimeter->overhang->... this condition is checked, to ensure the Arcs are generated within the ext perimeter only and not for each perimeter. Can be turned of to handle edge-case: only overhang-perimeter in gcode.
        "allowedSpaceForArcs": Polygon([[0,0],[500,0],[500,500],[0,500]]),#have control in which areas Arcs shall be generated
        "extendArcsIntoObject": 2, #extend the Arcs into the printed object for good adhesion. Unit:mm
        "plotStart":False, # plot the detected geoemtry in the prev Layer and the StartLine for Arc-Generation, use for debugging
        "plotArcs":False, #plot arcs for every filled polygon. use for debugging
        "CornerImportanceMultiplier":1, # Startpoint for Arc generation is chosen close to the middle of the StartLineString and at a corner. Higher=>Cornerselection more important.
        "minDistanceFromExtPerimeter":(gCodeSettingDict.get("perimeters")-1)*gCodeSettingDict.get("perimeter_extrusion_width"), #when having multiple perimeters the arcs would overlap with the inner ones. Regulate with this parameter,current=Toolpath can be coincident with perimeter,overlap=width/2
        "ArcPrintSpeed":3*60, #Unit:mm/min
        "ArcMinPrintSpeed":1*60,#Unit:mm/min
        "ArcSlowDownBelowThisDuration":3,# Arc Time below this Duration =>slow down, Unit: sec
        "ArcTravelFeedRate":50*60, # slower travel speed, Unit:mm/min
        "GCodeArcPtMinDist":0.1, # min Distance between points on the Arcs to for seperate GCode Command. Unit:mm
        "ArcCenterOffset":2, # Unit:mm, prevents very small Arcs by hiding the center in not printed section
        "ArcWidth":gCodeSettingDict.get("nozzle_diameter")*0.95, #change the spacing between the arcs,should be nozzle_diameter
        "ArcExtrusionMultiplier":1.35,
        "rMax":15, # the max radius of the arcs.
        "maxDistanceFromPerimeter":gCodeSettingDict.get("perimeter_extrusion_width")*1.5,#Control how much bumpiness you allow between arcs and perimeter. lower will follow perimeter better, but create a lot of very small arcs. Should be more that 1 Arcwidth! Unit:mm
        "pointsPerCircle":80, # each Arc starts as a discretized circle. Higher will slow down the code but give more accurate results for the arc-endings. 
        "extendArcDist":0.5 # extend Arcs for better bonding bewteen them, only end-piece affected(yet), Unit:mm
    }
    gCodeSettingDict.update(AddManualSettingsDict)
    return gCodeSettingDict

##### Helper Funcions to get Polygon from GCode #####
def getFileStreamAndPath(read=True):
    if len(sys.argv) != 2:
        print("Usage: python3 ex1.py <filename>")
        sys.exit(1)
    filepath = sys.argv[1]
    try:
        if read:
            f = open(filepath, "r")
        else:
            f=open(filepath, "w")    
        return f,filepath
    except IOError:
        print("File not found")
        sys.exit(1)
        
def splitGCodeIntoLayers(gcode:list)->list:
    gcode_list = []
    buff=[]
    for linenumber,line in enumerate(gcode):
        if ";LAYER_CHANGE" in line:
            gcode_list.append(buff)
            buff=[]
            buff.append(line)
        else:
            buff.append(line)
    gcode_list.append(buff)  #catch last layer
    print("last read linenumber:",linenumber)            
    return gcode_list
            

def getPtfromCmd(line:str)->Point:
    x=None
    y=None
    line=line.split(";")[0]
    cmds=line.split(" ")
    for c in cmds:
        if "X" in c:
            x=float(c[1:])
        elif "Y" in c:
            y=float(c[1:])
    if (x is not None) and (y is not None):
        p=Point(x,y)
    else:
        p=None    
    return p        

def makePolygonFromGCode(lines:list)->Polygon:
    pts=[]
    for line in lines:
        if ";WIPE" in line:
            break
        if "G1" in line:
            p=getPtfromCmd(line)
            if p:
                pts.append(p)
    if len(pts)>2:
        return Polygon(pts)
    else:
        #print("invalid poly: not enough pts")
        return None  

############ Classes #####        
class Layer():
    def __init__(self,lines=[],kwargs={})->None:
        self.lines=lines
        self.z=kwargs.get("z",None)
        self.polys=[]
        self.validpolys=[]
        self.extPerimeterPolys=[]
        self.binfills=[]
        self.features=[]
        self.oldpolys=[]
        self.dontPerformPerimeterCheck=kwargs.get('notPerformPerimeterCheck',False)
        self.deleteTheseInfills=[]
        self.deletelines=[]
        self.binfillsINFO=[]
    def extract_features(self)->None:
        buff=[]
        currenttype=""
        start=0
        for idl,line in enumerate(self.lines):
            if ";TYPE:" in line:
                if currenttype:
                    self.features.append([currenttype,buff,start])
                    buff=[]
                    start=idl                        
                currenttype=line
            else:
                buff.append(line)  
        self.features.append([currenttype,buff,start])# fetch last one

    def getRealFeatureStartPoint(self,idf:int)->str: # since GCode only stores destination of the move, the origin of the first move has to be included. 
        if idf<1:
            raise ValueError("Tried to find start point of first feature: not implemented")
        lines=self.features[idf-1][1]
        for line in reversed(lines):
            if "G1" in line:
                return line

    def makeOverhangPerimeter2Polys(self)->None:
        for idf,fe in enumerate(self.features):
            ftype=fe[0]
            lines=fe[1]
            if "Overhang" in ftype:
                if "External" in self.features[idf-1][0] or self.dontPerformPerimeterCheck:#overhang is in external perimeter. 
                    linesWithStart=[]
                    linesWithStart.append(self.getRealFeatureStartPoint(idf))
                    linesWithStart=linesWithStart+lines
                    poly=makePolygonFromGCode(linesWithStart)
                    if poly:
                        self.polys.append(poly)
    def makeExternalPerimeter2Polys(self)->None:
        extPerimeterIsStarted=False
        for idf,fe in enumerate(self.features):
            ftype=fe[0]
            lines=fe[1]
            if "External" in ftype or ("Overhang" in ftype and extPerimeterIsStarted) or ("Overhang" in ftype and self.dontPerformPerimeterCheck): #two different types of perimeter to for a poly: external perimeter and overhang perimeter + option for manual errorhandling, when there is no feature "external"
                if not extPerimeterIsStarted:
                    linesWithStart=[]
                    linesWithStart.append(self.getRealFeatureStartPoint(idf))
                linesWithStart=linesWithStart+lines
                extPerimeterIsStarted=True
            elif extPerimeterIsStarted or (idf==len(self.features)-1 and extPerimeterIsStarted):#finish the poly if other feature or end of featurelist
                poly=makePolygonFromGCode(linesWithStart)
                if poly:
                    self.extPerimeterPolys.append(poly) 
                extPerimeterIsStarted=False      
    def makeStartLineString(self,startpt:Point,endpt:Point,kwargs={})->LineString:
        if not self.extPerimeterPolys:
            self.makeExternalPerimeter2Polys()
        startOnPoly=None
        endOnPoly=None
        startDist=1e4
        endDist=1e4
        polyIndex=0
        for idp,poly in enumerate(self.extPerimeterPolys):
            start=nearest_points(poly,startpt)[0]
            end=nearest_points(poly,endpt)[0]
            if startDist>start.distance(startpt):
                startDist=start.distance(startpt)
                polyIndex=idp
        poly=self.extPerimeterPolys[polyIndex]
        if kwargs.get("plotStart"): 
            plot_geometry(poly)
        poly=poly.buffer(-1*kwargs.get("extendArcsIntoObject",2)) #erode for overlap
        startOnPoly=nearest_points(poly,startpt)[0]   #nearest_points is enough, because overhang always starts at some kind of a corner.  
        endOnPoly=nearest_points(poly,endpt)[0]
        lineHasStarted=False
        startLine=[]
        ids=0
        pts=[Point(p) for p in poly.exterior.coords]
        while ids<len(poly.exterior.coords)*2:
            p=pts[ids%len(pts)]
            #print("point:",p,"StartDist:",startOnPoly.distance(p))
            if startOnPoly.distance(p)<1e-2:
                lineHasStarted=True
            if lineHasStarted:
                startLine.append(p)
            if endOnPoly.distance(p)<1e-2 and lineHasStarted:
                break
            ids+=1
        if len(startLine)<1:
            raise ValueError("StartLine not found")    
        ls=LineString(startLine)
        if kwargs.get("plotStart"):
            #plot_geometry([Point(p) for p in poly.exterior.coords],changecolor=True)
            plot_geometry(poly,color="b")
            #plot_geometry(startOnPoly,color="g")
            #plot_geometry(endOnPoly,color="y")
            plot_geometry(ls,color="m")
            plt.title("Start-Geometry in prev Layer")
            plt.legend(["ExtPerimeter","Eroded Perimeter","StartLine for Arc Generation"])
            plt.show()
        return ls

    def spotBridgeInfill(self)->None:
        for idf,fe in enumerate(self.features):
            ftype=fe[0]
            lines=fe[1]
            start=fe[2]        
            if "Bridge infill" in ftype:
                infillpts=[]
                for line in lines:
                    if "G1" in line:
                        p=getPtfromCmd(line)
                        if p:
                            infillpts.append(p)
                self.binfills.append(infillpts)
                self.binfillsINFO.append(idf)   
    
    def verifypolys(self)->None: #verify poly is inside allowed space and having at least one pt of bridging infill inside it.
        for poly in self.polys:
            valid=False
            allowedSpacePolygon=parameters.get("allowedSpaceForArcs")
            if allowedSpacePolygon:
                if allowedSpacePolygon.contains(poly):
                    for idb,binfill in enumerate(self.binfills):
                        for p in binfill:
                            if poly.contains(p):
                                valid=True
                                self.validpolys.append(poly)
                                self.deleteTheseInfills.append(idb)
                                break
                        if valid:
                            break
    def prepareDeletion(self)->None:
        if self.deleteTheseInfills:
            for idd in self.deleteTheseInfills:
                idf=self.binfillsINFO[idd]
                start=self.features[idf][2]
                if idf<len(self.features)-1:
                    end=self.features[idf+1][2]
                else:
                    end=len(self.lines)    
                self.deletelines.append([start,end])
    def exportThisLine(self,linenumber:int)->bool:
        export=True
        if len(self.deletelines)>0:
            for d in self.deletelines:
                if linenumber>=d[0] and linenumber<=d[1]:
                    export=False
        return export            

                            
class Arc():
    def __init__(self,center:Point,r:float,kwargs={}) -> None:
        self.center=center
        self.r=r
        self.pointsPerCircle=kwargs.get("pointsPerCircle",80)
    def setPoly(self,poly:Polygon)->None:
        self.poly=poly    
    def extractArcBoundary(self)->LineString:
        circ=create_circle(self.center,self.r,self.pointsPerCircle)    
        trueArc=self.poly.boundary.intersection(circ.boundary.buffer(1e-2))
        if trueArc.geom_type=='MultiLineString':
            merged=linemerge(trueArc)
        elif trueArc.geom_type=='LineString':
            self.arcline=trueArc
            return trueArc
        else:
            #print("Other Geom-Type:",trueArc.geom_type)
            merged=linemerge(MultiLineString([l for l in trueArc.geoms if l.geom_type=='LineString']))
        if merged.geom_type=="LineString":
            self.arcline=merged
            return merged                             
        elif merged.geom_type=="MultiLineString": #Error handling
            warnings.warn("Boundary was not continous, only first part passed")
            self.arcline=merged.geoms[0]
            return merged.geoms[0]
        else:
            raise ValueError("ArcBoundary merging Error")
    def generateConcentricArc(self,startpt:Point,remainingSpace:Polygon|MultiPolygon)->Polygon:
        circ=create_circle(startpt,self.r,self.pointsPerCircle)
        arc=circ.intersection(remainingSpace)
        self.poly=arc
        return arc            

############ Helper Functions Arc Generation #####
def midpoint(p1:Point, p2:Point):
    return Point((p1.x + p2.x)/2, (p1.y + p2.y)/2)

def getstartptOnLS(ls:LineString,kwargs={})->Point:
    if len(ls.coords)<2:
        raise ValueError("Start LineString with <2 Points invalid")
    if len(ls.coords)==2:
        return midpoint(Point(ls.coords[0]),Point(ls.coords[1]))
    scores=[]
    curLength=0
    pts=[Point(p) for p in ls.coords]
    coords=[np.array(p) for p in ls.coords]
    for idp,p in enumerate(pts):
        if idp==0 or idp==len(pts)-1:
            scores.append(0)
            continue
        curLength+=p.distance(pts[idp-1])
        relLength=curLength/ls.length
        lengthscore=1-np.abs(relLength-0.5)#Hat-function: pointscore=1 at relLength=0.5, 0 at start or end.
        v1=coords[idp]-coords[idp-1]
        v2=coords[idp+1]-coords[idp]
        if np.linalg.norm(v1)>0 and np.linalg.norm(v2)>0:#calc angle only for non-zero-vectors
            v1=v1/np.linalg.norm(v1)
            v2=v2/np.linalg.norm(v2)
            anglescore=np.abs(np.sin(np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))))#prefer points at corners
            anglescore*=kwargs.get("CornerImportanceMultiplier",1)
            scores.append(lengthscore+anglescore)
        else:
            scores.append(lengthscore)    
    maxIndex=scores.index(max(scores))
    return pts[maxIndex]

def create_circle(p:Point, radius:float, n:int)->Polygon:
    x=p.x
    y=p.y
    return Polygon([[radius*np.sin(theta)+x, radius*np.cos(theta)+y] for theta in np.linspace(0, 2*np.pi - 2*np.pi/n, int(n))])       

def get_farthest_point(arc:Polygon, base_poly:Polygon, remaining_empty_space:Polygon):#function ported from Steven McCulloch
    """
    Find the point on a given arc that is farthest away from the base polygon.
    In other words, the point on which the largest circle can be drawn without going outside the base polygon.
    Parameters
    ----------
    arc: Polygon
        The arc in question
    base_poly: Polygon
        The base polygon
    remaining_empty_space: Polygon
        The polygon representing the space left to be filled in the base polygon        
    Returns
    -------
    farthest_point: Point
        The point on the arc that is farthest away from the base polygon
    longest_distance: float
        How far away the polygon is from the farthest point
    point_on_poly: Point
        The point on the base polygon that is closest to the arc
    """
    longest_distance = -1
    farthest_point = Point([0, 0])
    pointFound=False
    # Handle input for polygons and LineString
    # The first arc begins on a LineString rather than a Polygon
    if arc.geom_type == 'Polygon':
        arc_coords = arc.exterior.coords
    elif arc.geom_type == 'LineString':
        arc_coords = np.linspace(list(arc.coords)[0], list(arc.coords)[1])
    else:
        print('get_farthest_distance: Wrong shape type given',type(arc))
        plt.title("Function get_farthest_point went wrong")
        plot_geometry(base_poly,"b")
        plot_geometry(arc,"r")
        plt.show()
    # For every point in the arc, find out which point is farthest away from the base polygon
    for p in list(arc_coords):
        distance = Point(p).distance(base_poly.boundary)
        if (distance > longest_distance) and ((remaining_empty_space.buffer(1e-2).contains(Point(p)))):
            longest_distance = distance
            farthest_point = Point(p)
            pointFound = True
    point_on_poly = nearest_points(base_poly, farthest_point)[0]
    if pointFound:
        return farthest_point, longest_distance, point_on_poly 
    else:
        return None, None, None       

def move_toward_point(start_point:Point, target_point:Point, distance:float)->Point:
    """Moves a point a set distance toward another point"""

    # Calculate the direction in which to move
    dx = target_point.x - start_point.x
    dy = target_point.y - start_point.y

    # Normalize the direction
    magnitude = (dx**2 + dy**2)**0.5
    dx /= magnitude
    dy /= magnitude

    # Move the point in the direction of the target by the set distance
    return Point(start_point.x + dx*distance, start_point.y + dy*distance)        


    


    
def generateMultipleConcentricArcs(startpt:Point,rMin:float,rMax:float, boundaryLineString:LineString,remainingSpace:Polygon|MultiPolygon,kwargs={})->list:
    arcs=[]
    r=rMin
    while r<=rMax:
        arcObj=Arc(startpt,r,kwargs=kwargs)
        arc=arcObj.generateConcentricArc(startpt,remainingSpace)
        if arc.intersects(boundaryLineString):
            break
        arcs.append(arcObj)
        #print("True Arc type:",type(arc4gcode))
        r+=kwargs.get("ArcWidth")
    return arcs


############ Helper Functions Arc Validation #####

def getValueBasedColor(val:float, max_val=10)->tuple:
    normalizedVal = val / max_val
    rgb=[0,0,0]
    rgb[0]=min(normalizedVal,1)
    rgb[2]=1-rgb[0]
    return tuple(rgb)


def plot_geometry(geometry, color='black', linewidth=1,**kwargs):
    if type(geometry)==type([]):
        for idx,geo in enumerate(geometry):
            if kwargs.get("changecolor"):
                color=getValueBasedColor(idx,len(geometry))
            plot_geometry(geo,color=color,linewidth=linewidth,kwargs=kwargs)
    elif geometry.geom_type == 'Point':
        x, y = geometry.x, geometry.y
        plt.scatter(x, y, color=color, linewidth=linewidth)        
    elif geometry.geom_type == 'LineString':
        x, y = geometry.xy
        plt.plot(x, y, color=color, linewidth=linewidth)
    elif geometry.geom_type == 'Polygon':
        x, y = geometry.exterior.xy
        plt.plot(x, y, color=color, linewidth=linewidth)
        if kwargs.get("filled"):
            plt.fill(x,y,color=color,alpha=0.8)
        for interior in geometry.interiors:
            x, y = interior.xy
            plt.plot(x, y, color=color, linewidth=linewidth)
            if kwargs.get("filled_holes"):
                plt.fill(x,y,color=color,alpha=0.5)
    elif geometry.geom_type == 'MultiLineString':
        for line in geometry.geoms:
            x, y = line.xy
            plt.plot(x, y, color=color, linewidth=linewidth)
    elif geometry.geom_type == 'MultiPolygon':
        for polygon in geometry.geoms:
            plot_geometry(polygon,color=color,linewidth=linewidth,kwargs=kwargs)
    else:
        print('Unhandled geometry type: ' + geometry.geom_type)

        

############ Helper Functions write finished GCode #####
def readSettingsFromGCode2dict(gcodeLines:list)->dict:
    gCodeSettingDict={}
    isSetting=False
    for line in gcodeLines:
        if "; prusaslicer_config = begin" in line:
            isSetting=True
            continue
        if isSetting :
            setting=line.strip(";").strip("\n").split("= ")
            if len(setting)==2:
                try:
                    gCodeSettingDict[setting[0].strip(" ")]=literal_eval(setting[1]) # automaticly convert into int,float,...
                except:
                    gCodeSettingDict[setting[0].strip(" ")]=setting[1] # leave the complex settings as strings. They shall be handled individually if necessary 
            else:
                print("unmatching setting:",setting)          
    return gCodeSettingDict
def checkforNecesarrySettings(gCodeSettingDict:dict)->bool:
    if not ";AFTER_LAYER_CHANGE" in gCodeSettingDict.get("layer_gcode"):
        warnings.warn("After Layer Change missing keyword, expected: ';AFTER_LAYER_CHANGE' ")
        return False
    if not gCodeSettingDict.get("use_relative_e_distances"):
        warnings.warn("Script only works with relative e-distances. Change acordingly.")
        return False
    if gCodeSettingDict.get("infill_first"):
        warnings.warn("Infill is printed before perimeter. This can cause problems with the script.")
    if gCodeSettingDict.get("external_perimeters_first"):
        warnings.warn("External perimeter is printed before inner perimeters. Change for better overhang performance. ")    
    return True
def calcEStepsPerMM(settingsdict:dict)->float:
    eVol = (settingsdict.get("nozzle_diameter")/2)**2 * np.pi *settingsdict.get("ArcExtrusionMultiplier",1)#printing in midair will result in circular shape. Scource: https://manual.slic3r.org/advanced/flow-math
    if settingsdict.get("use_volumetric_e"):
        return eVol
    else:
        eInMm=eVol/((settingsdict.get("filament_diameter")/2)**2 *np.pi)
        return eInMm


def p2GCode(p:Point,E=0,**kwargs)->str:
    line=f"G1 X{p.x:.4} Y{p.y:.4} "
    line+="E0" if E==0 else f"E{E:.5f}"
    if kwargs.get('F'):
        line+=f" F{kwargs.get('F'):0d}"
    line+='\n'       
    return line  

def retractGCode(retract=True,kwargs={})->str:
    retractDist=kwargs.get("retract_length",1)
    E= -retractDist if retract else retractDist
    return f"G1 E{E} F{kwargs.get('retract_speed',2100)}\n"  

def setFeedRateGCode(F:int)->str:
    return f"G1 F{F}\n"     

def arc2GCode(arcline:LineString,eStepsPerMM:float,arcidx=None,kwargs={})->list:
    GCodeLines=[]
    p1=None
    pts=[Point(p) for p in arcline.coords]
    #plt.plot([p.x for p in pts],[p.y for p in pts])
    #plt.show()      
    extDist=kwargs.get("extendArcDist",0.5)
    pExtend=move_toward_point(pts[-2],pts[-1],extDist)
    arcPrintSpeed=np.clip(arcline.length/(kwargs.get("ArcSlowDownBelowThisDuration",3))*60,
                            kwargs.get("ArcMinPrintSpeed",1*60),kwargs.get('ArcPrintSpeed',2*60)) # *60 bc unit conversion:mm/s=>mm/min
    for idp,p in enumerate(pts):
        if idp==0:
            p1=p
            GCodeLines.append(f";Arc {arcidx if arcidx else ' '} Length:{arcline.length}\n")
            GCodeLines.append(p2GCode(p,F=kwargs.get('ArcTravelFeedRate',100*60)))#feedrate is mm/min...
            GCodeLines.append(retractGCode(retract=False))
            GCodeLines.append(setFeedRateGCode(arcPrintSpeed))
        else:
            dist=p.distance(p1)
            if dist>kwargs.get("GCodeArcPtMinDist",0.1):
                GCodeLines.append(p2GCode(p,E=dist*eStepsPerMM))
                p1=p
        if idp==len(pts)-1:
            GCodeLines.append(p2GCode(pExtend,E=extDist*eStepsPerMM))#extend arc tangentially for better bonding between arcs
            GCodeLines.append(retractGCode(retract=True))
    return GCodeLines        


############# MAIN ############
if __name__=="__main__":
    gCodeFileStream,path2GCode = getFileStreamAndPath()
    gCodeLines=gCodeFileStream.readlines()
    gCodeSettingDict=readSettingsFromGCode2dict(gCodeLines)
    parameters=makeFullSettingDict(gCodeSettingDict)
    if not checkforNecesarrySettings(gCodeSettingDict):
        raise ValueError("Incompatible Settings used!")
    layerobjs=[]
    if gCodeFileStream:
        layers=splitGCodeIntoLayers(gCodeLines)
        gCodeFileStream.close()
        print("layers:",len(layers))
        for idl,layerlines in enumerate(layers):
            layer=Layer(layerlines,parameters)
            layerobjs.append(layer)   
        for idl,layer in enumerate(layerobjs):    
            if idl<1:
                continue # no overhangs in the first layer and dont mess with the setup
            else:
                layer.extract_features()
                layer.makeOverhangPerimeter2Polys()
                if layer.polys:
                    print(f"overhang found layer {idl}:",len(layer.polys))
                    layer.spotBridgeInfill()
                    layer.verifypolys()
                    #print("validpolys:",len(layer.validpolys))
                    for poly in layer.validpolys:
                        #make parameters more readable
                        maxDistanceFromPerimeter=parameters.get("maxDistanceFromPerimeter") # how much 'bumpiness' you accept in the outline. Lower will generate more small arcs to follow the perimeter better (corners!). Good practice: 2 perimeters+ threshold of 2width=minimal exact touching (if rMin satisfied)
                        rMax=parameters.get("rMax",15)
                        pointsPerCircle=parameters.get("pointsPerCircle",80)
                        arcWidth=parameters.get("ArcWidth")
                        rMin=parameters.get("ArcCenterOffset")+arcWidth/1.5
                        #initialize
                        finalarcs=[]
                        arcs=[]
                        arcs4gcode=[]
                        pts=[Point(p) for p in poly.exterior.coords]
                        #find StartPoint and StartLineString
                        prevLayer=layerobjs[idl-1]
                        prevLayer.makeExternalPerimeter2Polys()
                        startLineString=prevLayer.makeStartLineString(pts[-2],pts[-1],parameters)
                        boundaryLineString=LineString(poly.exterior.coords[:-1])
                        thresholdedboundaryLS=boundaryLineString.buffer(parameters.get("minDistanceFromExtPerimeter",0.25))
                        polyWithCompleteBoundary=Polygon(list(startLineString.coords)+list(boundaryLineString.coords))
                        startpt=getstartptOnLS(startLineString,parameters)
                        remainingSpace=polyWithCompleteBoundary
                        #first step in Arc Generation
                        concentricArcs=generateMultipleConcentricArcs(startpt,rMin,rMax,thresholdedboundaryLS,remainingSpace,parameters)
                        arcBoundarys=[arc.extractArcBoundary() for arc in concentricArcs]
                        #print(f"number of concentric arcs generated:",len(concentricArcs))
                        if concentricArcs:
                            finalarcs.append(concentricArcs[-1]) 
                            for arc in concentricArcs: 
                                remainingSpace=remainingSpace.difference(arc.poly.buffer(1e-2))
                                arcs.append(arc)
                            for arcboundary in arcBoundarys:    
                                arcs4gcode.append(arcboundary)
                        else:
                            raise ValueError("Initialization Error")
                       
                        #start bfs (breadth first search algorithm) to fill the remainingspace
                        idx=0
                        safetyBreak=0
                        while idx<len(finalarcs):
                            print("while executed:",idx, len(finalarcs))
                            curArc=finalarcs[idx]
                            farthestPointOnArc,longestDistance,NearestPointOnPoly=get_farthest_point(curArc.poly,polyWithCompleteBoundary,remainingSpace)
                            if not farthestPointOnArc or longestDistance<maxDistanceFromPerimeter:#no more pts on arc
                                idx+=1 #go to next arc
                                continue
                            startpt=move_toward_point(farthestPointOnArc,curArc.center,parameters.get("ArcCenterOffset",2))
                            concentricArcs=generateMultipleConcentricArcs(startpt,rMin,rMax,thresholdedboundaryLS,remainingSpace,parameters)
                            arcBoundarys=[arc.extractArcBoundary() for arc in concentricArcs]
                            #print(f"number of concentric arcs generated:",len(concentricArcs))
                            if len(concentricArcs)>0:
                                for arc in concentricArcs: 
                                    remainingSpace=remainingSpace.difference(arc.poly.buffer(1e-2))
                                    arcs.append(arc)
                                finalarcs.append(concentricArcs[-1])
                                for arcboundary in arcBoundarys:    
                                    arcs4gcode.append(arcboundary)
                            else:
                                idx+=1 # no possible concentric arcs found= arc complete, proceed to next
                            safetyBreak+=1
                            if safetyBreak>2000:
                                break
                            if parameters.get("plotArcs"):
                                plt.title(f"Iteration {idx}, Total No Start Points: {len(finalarcs)}, Total No Arcs: {len(arcs)}")
                                plot_geometry(startLineString,'r')
                                plot_geometry(boundaryLineString)
                                plot_geometry([arc.poly for arc in arcs],changecolor=True)
                                plot_geometry(remainingSpace,'g',filled=True)
                                plot_geometry(startpt,"r")
                                plt.axis('square')
                                plt.show()  
                        
                        #generate gcode for arc and insert at the beginning of the layer
                        eStepsPerMM=calcEStepsPerMM(parameters)
                        arcOverhangGCode=[]
                        arcOverhangGCode.append(f"M106 S{np.round(parameters.get('bridge_fan_speed',100)*2.55)}")#turn cooling Fan on at Bridge Setting
                        #for arc in arcs4gcode:
                        #    plot_geometry(arc)
                        #    plot_geometry(Point(arc.coords[0]))
                        #plt.show()
                        for ida,arc in enumerate(arcs4gcode):    
                            arcGCode=arc2GCode(arcline=arc,eStepsPerMM=eStepsPerMM,arcidx=ida)
                            arcOverhangGCode.append(arcGCode)
                        modifiedlayer=Layer() # copy the other infos if needed: future to do
                        isInjected=False
                        layer.prepareDeletion()
                        #print("FEATURES:",[(f[0],f[2]) for f in layer.features])
                        injectionStart=None
                        for idline,line in enumerate(layer.lines):
                            if ";TYPE" in line and not isInjected:
                                injectionStart=idline
                                modifiedlayer.lines.append(";TYPE:Bridge infill\n")
                                for overhangline in arcOverhangGCode:
                                    for arcline in overhangline:
                                        for cmdline in arcline:
                                            modifiedlayer.lines.append(cmdline)
                                isInjected=True
                                #add restored pre-injected tool position
                                for id in reversed(range(injectionStart)):
                                    if "X" in layer.lines[id]:
                                        modifiedlayer.lines.append(layer.lines[id])
                                        break
                            if layer.exportThisLine(idline):        
                                modifiedlayer.lines.append(line)
                        layerobjs[idl]=modifiedlayer  # overwrite the infos
f=open(path2GCode,"w")
print("write to file")
for layer in layerobjs:
    f.writelines(layer.lines)
f.close()    
#os.startfile(path2GCode, 'open')
print("Code executed successful")
input("Push enter to close this window")
