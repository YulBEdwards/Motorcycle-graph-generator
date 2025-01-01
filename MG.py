import matplotlib.pyplot as plt
import math

# Line segment class
class Segment:
    def __init__(self, sid, p1, p2, a, b, c):
        self.sid = sid  # Unique identifier
        self.p1  = p1   # (x1, y1)
        self.p2  = p2   # (x2, y2)
        self.a   = a    # line is a*x+b*y+c
        self.b   = b
        self.c   = c

class LineDef:
    def __init__(self, a, b, c):
        self.a   = a    # line is a*x+b*y+c
        self.b   = b
        self.c   = c

class Track:
    def __init__(self, tid, sid, t, speed, bisect, target):
        self.tid      = tid     # Track identifier
        self.sid      = sid     # Corner identifier
        self.t        = t       # Time at ending intersection
        self.speed    = speed
        self.bisect   = bisect  # Bisector definition for track
        self.target   = target  # Side id of target        
    
class Corner:
    def __init__(self, cid, bisect):
        self.cid    = cid
        self.bisect = bisect
        
# Track intersect object
class TI:
    def __init__(self, idsml, idlrg, tsml, tlrg, xi, yi):
        self.idsml = idsml
        self.idlrg = idlrg
        self.tsml  = tsml
        self.tlrg  = tlrg
        self.xi    = xi
        self.yi    = yi
        
def generateSegments():

    s = [[666, -134],[603, -61], [491, -40], [447, -95],
         [497, -160],[420, -179],[300, -105],[255, -123],
         [157, -68], [101, -158],[149, -217],[108, -323],
         [226, -379],[201, -470],[238, -542],[315, -525],
         [407, -553],[450, -441],[527, -472],[531, -364],
         [629, -365],[599, -252],[675, -213],[666, -134]]
   
    Segments = []         
    numSegments = len(s)-1
    clockwise = 0
    for i in range(numSegments+2):
        if i>numSegments-1:
            ix = i-numSegments
        else:
            ix=i
        p1 = s[ix]
        p2 = s[ix+1]
        if i<numSegments:
            cross = p1[0]*p2[1]-p1[1]*p2[0]
            clockwise = clockwise+cross
        a  = p2[1]-p1[1]
        b  = p1[0]-p2[0]
        xm = 0.5*(p1[0]+p2[0])
        ym = 0.5*(p1[1]+p2[1])
        c  = -a*xm-b*ym
        L  = math.sqrt(a**2+b**2)
        if L>0:
            a = a/L
            b = b/L
            c = c/L
        Segments.append(Segment(sid=i, p1=p1, p2=p2, a=a, b=b, c=c))
        
    if clockwise<0:
        print("closes clockwise")
    else:
        print("closes counter clockwise")
        
    return Segments, numSegments

def generateCorners(Segments, numSegments):
    corners = []
    tracks  = []
    num_Corners = numSegments   
    clockwise = 0
    for i in range(num_Corners):
        s1 = Segments[i]
        s2 = Segments[i+1]
        turning = s1.a*s2.b-s1.b*s2.a
        reflex = turning<0
        # determine bisector
        a = s2.a-s1.a
        b = s2.b-s1.b
        c = s2.c-s1.c
        L  = math.sqrt(a**2+b**2)
        if L>0:
            a = a/L
            b = b/L
            c = c/L
            bisect = Segment(sid=i, p1=s1.p2, p2=[0,0], a=a, b=b, c=c)   
        if not reflex:          
            Segments[i].cid = len(corners)
            Segments[i].tid = -1
            corners.append(Corner(cid=i, bisect=bisect))
        else:
            tid = len(tracks)
            Segments[i].tid = len(tracks)
            Segments[i].cid = -1        
            cos_turning = s1.a*s2.a+s1.b*s2.b
            speed = math.sin(0.5*(math.acos(cos_turning)+math.pi))            
            tracks.append(Track(tid = tid, sid=i,  \
                t = 0, speed=speed, bisect=bisect, target=0))
    return corners, tracks, Segments

# find target of a ray starting at vertex corner
def intersectionSegmentTrack(segment,numSegments,track):
    if track.sid==numSegments-1:
        ix = 0
    else:
        ix = track.sid+1
    if segment.sid == track.sid or segment.sid == ix:        
        return False, None, None, None
    x1 = segment.p1[0]
    y1 = segment.p1[1]
    x2 = segment.p2[0]
    y2 = segment.p2[1]
    x0 = track.bisect.p1[0]
    y0 = track.bisect.p1[1]
    dx = track.bisect.b
    dy = -track.bisect.a
    delta = dx*(y2-y1)-dy*(x2-x1)
    t = ((x1-x0)*(y2-y1)-(y1-y0)*(x2-x1))/delta
    u = -(dx*(y1-y0)-dy*(x1-x0))/delta
    intersect  = u>=0 and u<=1 and t>=0
    xi = x0+t*dx
    yi = y0+t*dy
    target = segment.sid    
    # get the closest segment to the base point
    if intersect and track.bisect.p2 != [0,0]:
        d1 = (track.bisect.p1[0]-track.bisect.p2[0])**2+ \
             (track.bisect.p1[1]-track.bisect.p2[1])**2
        d2 = (track.bisect.p1[0]-xi)**2+ \
             (track.bisect.p1[1]-yi)**2
        if d1<d2:         
            return False, None, None, None
    return intersect,xi,yi,target

def intersectionTrackTrack(track1,track2):
    x1 = track1.bisect.p1[0]
    y1 = track1.bisect.p1[1]
    x2 = track1.bisect.p2[0]
    y2 = track1.bisect.p2[1]
    x3 = track2.bisect.p1[0]
    y3 = track2.bisect.p1[1]
    x4 = track2.bisect.p2[0]
    y4 = track2.bisect.p2[1]
    t = ((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/  \
        ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
    u = -((x1-x2)*(y1-y3)-(y1-y2)*(x1-x3))/  \
        ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
    intersect  = u<=1 and u>=0 and t<=1 and t>=0
    xi = x1+t*(x2-x1)
    yi = y1+t*(y2-y1)
    return intersect,xi,yi

def distanceToEdge(p,segment):
    return abs(p[0]*segment.a+p[1]*segment.b+segment.c)

def intersectionLineLine(line1,line2):
    delta = lin1.a*line2.b - line2.a*line1.b
    x = (line1.b*line2.c-line2.b*line1.c)/delta
    y = (line2.a*line1.c-line1.a*line2.c)/delta

def generateBisector(seg1,seg2):
    a = seg2.a-seg1.a
    b = seg2.b-seg1.b
    c = seg2.c-seg1.c
    L  = math.sqrt(a**2+b**2)
    if L>0:
        a = a/L
        b = b/L
        c = c/L
    return LineDef(a=a, b=b, c=c)           


def intersectionBisectBisect(seg1A,seg1B,seg2A,seg2B):
    LD1 = generateBisector(Segments[seg1A],Segments[seg1B])
    LD2 = generateBisector(Segments[seg2A],Segments[seg2B])
    delta = LD1.a*LD2.b-LD2.a*LD1.b
    x = (LD1.b*LD2.c-LD2.b*LD1.c)/delta
    y = (LD2.a*LD1.c-LD1.a*LD2.c)/delta
    return x,y


# Set up the bounds
bounds = (0, -600, 800, 0)

# Insert line Segments
Segments, numSegments = generateSegments()

# determine the corner bisectors and tracks
Corners, Tracks, Segments = generateCorners(Segments, numSegments)
    
# Setup the graph
fig, ax = plt.subplots(figsize=(8, 8))

# Identify the corners
for segment in Segments:
    x1, y1 = segment.p1
    x2, y2 = segment.p2
    ax.plot([x1, x2], [y1, y2], 'r-')
    if segment.sid<numSegments:
        ax.text(0.5*(x1+x2),0.5*(y1+y2), str(segment.sid), \
         color="black", fontsize=7)
    #ax.text(x1,y1, str(segment.sid), \
    #    color="black", fontsize=7)

# Determine intersection of reflex ray and closest corner
for track in Tracks:
    for i in range(numSegments):
        segment = Segments[i]
        cross, xi, yi, target =  \
         intersectionSegmentTrack(segment, numSegments, track)
        if cross:
            track.bisect.p2[0]=xi
            track.bisect.p2[1]=yi            
            track.target = target
    
# Generate the intersection time pairs        
TIs = []

for track1 in Tracks:
    for track2 in Tracks:
        if track1.sid > track2.sid:
            cross,xi,yi = intersectionTrackTrack(track1,track2)            
            if cross:
                d1 = math.sqrt((xi-track1.bisect.p1[0])**2 + \
                           (yi-track1.bisect.p1[1])**2)
                d2 = math.sqrt((xi-track2.bisect.p1[0])**2 + \
                           (yi-track2.bisect.p1[1])**2)
                t1 = d1 * track1.speed
                t2 = d2 * track2.speed
                i = len(TIs)+1
                if t1<=t2:
                    TIs.append(TI(idsml=track1.tid, idlrg=track2.tid,\
                        tsml=t1, tlrg=t2, xi=xi, yi=yi))
                else:
                    TIs.append(TI(idsml=track2.tid, idlrg=track1.tid, \
                        tsml=t2, tlrg=t1, xi=xi, yi=yi))
                    
print("Number of edges:", numSegments)
L = len(Tracks)                    
print("Candidate pairs:", int(L*(L+1)/2))
print("Pairs used:",len(TIs))
      
# Sort the time pairs                    
def get_tsml(x):
    return x.tsml
            
TIs.sort(key=get_tsml)

# Determine death location
for TI in TIs:
    if(TI.tlrg < Tracks[TI.idlrg].t and (TI.tsml<Tracks[TI.idsml].t or \
            Tracks[TI.idsml].t==0) or Tracks[TI.idlrg].t==0):
        cross,xi,yi =  cross,xi,yi = \
            intersectionTrackTrack(Tracks[TI.idsml],Tracks[TI.idlrg])
        if cross:
            Tracks[TI.idlrg].t = TI.tlrg
            Tracks[TI.idlrg].bisect.p2[0] = TI.xi
            Tracks[TI.idlrg].bisect.p2[1] = TI.yi
    
# plot the tracks
for track in Tracks:
    x1, y1 = track.bisect.p1
    x2, y2 = track.bisect.p2
    ax.plot([x1, x2], [y1, y2], color = "blue", linewidth=0.75 )
    ax.plot([x2],[y2], 'o', color ="blue", ms=1.75)
    
ax.set_aspect('equal', 'box')
ax.set_xlim(bounds[0], bounds[2])
ax.set_ylim(bounds[1], bounds[3])
plt.title("Motorcycle graph")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()
     
