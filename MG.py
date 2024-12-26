import matplotlib.pyplot as plt
import math

# Line segment class
class Segment:
    def __init__(self, sid, p1, p2, a, b, c):
        self.sid = sid # Unique identifier
        self.p1  = p1  # (x1, y1)
        self.p2  = p2  # (x2, y2)
        self.a   = a   # line is a*x+b*y+c
        self.b   = b
        self.c   = c

class Track:
    def __init__(self, tid, sid, t, speed, bisect):
        self.tid      = tid     # Track identifier
        self.sid      = sid     # Bisect identifier
        self.t        = t       # Time at ending intersection
        self.speed    = speed
        self.bisect   = bisect  # Bisector definition for track
    
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
    s = [[666, -134],[603, -61], [491, -40], [447, -95],  \
         [487, -167],[398, -179],[339, -105],[255, -123], \
         [157, -68], [101, -138],[149, -227],[108, -323], \
         [226, -379],[201, -470],[238, -542],[315, -525], \
         [427, -553],[439, -451],[527, -472],[531, -364], \
         [639, -355],[599, -252],[675, -213],[666, -134]]
    Segments = []
    num_Segments = len(s)-1
    
    for i in range(num_Segments):
            p1 = s[i]
            p2 = s[i+1]            
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
            Segments.append(Segment(sid=i+1, p1=p1, p2=p2, \
               a=a, b=b, c=c))                
    return Segments

def generateCorners(Segments):
    corners = []
    tracks  = []
    num_Segments = len(Segments)-1
    for i in range(num_Segments):
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
            bisect = Segment(sid=i+1, p1=s1.p2, p2=[0,0], a=a, b=b, c=c)   
        if not reflex:          
            corners.append(Corner(cid=i+1, bisect=bisect))
        else:
            tid = len(tracks)+1
            cos_turning = s1.a*s2.a+s1.b*s2.b
            speed = math.sin(0.5*(math.acos(cos_turning)+math.pi))            
            tracks.append(Track(tid = tid, sid=i+1,  \
                t = 0, speed=speed, bisect=bisect))            
    return corners, tracks
"""
def intersectionTrackTrack(track1,track2):
    Aa = track1.bisect.a
    Ab = track1.bisect.b
    Ac = track1.bisect.c
    Ba = track2.bisect.a
    Bb = track2.bisect.b
    Bc = track2.bisect.c
    delta = Aa*Bb-Ba*Ab
    x = (Ab*Bc-Bb*Ac)/delta
    y = (Ac*Ba-Bc*Aa)/delta    
    return x,y
"""
def intersectionSegmentTrack(segment,track):
    if segment.sid == track.sid or segment.sid == track.sid+1:
        return False, None, None
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
    # get the closest segment to the base point
    if intersect and track.bisect.p2 != [0,0]:
        d1 = (track.bisect.p1[0]-track.bisect.p2[0])**2+ \
             (track.bisect.p1[1]-track.bisect.p2[1])**2
        d2 = (track.bisect.p1[0]-xi)**2+ \
             (track.bisect.p1[1]-yi)**2
        if d1<d2:         
            return False, None, None
    return intersect,xi,yi

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

# Set up the bounds
bounds = (0, -600, 800, 0)

# Insert random line Segments
Segments = generateSegments()

# determine the corner bisectors and tracks
Corners, Tracks = generateCorners(Segments)

# Setup the graph
fig, ax = plt.subplots(figsize=(8, 8))

# Identify the corners
for segment in Segments:
    x1, y1 = segment.p1
    x2, y2 = segment.p2
    ax.plot([x1, x2], [y1, y2], 'r-')
    ax.text(x2,y2, str(segment.sid), \
         color="black", fontsize=6)

# Determine intersection of reflex ray and closest corner
for track in Tracks:
    for segment in Segments:
        cross, xi, yi = intersectionSegmentTrack(segment, track)
        if cross:
            track.bisect.p2[0]=xi
            track.bisect.p2[1]=yi            
        
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
                    
print("Number of edges:", len(Segments))
L = len(Tracks)                    
print("Candidate pairs:", int(L*(L+1)/2))
print("Pairs used:",len(TIs))
      
# Sort the time pairs                    
def get_tsml(x):
    return x.tsml
            
TIs.sort(key=get_tsml)

# Determine death location
for TI in TIs:
    if(TI.tlrg < Tracks[TI.idlrg-1].t and (TI.tsml<Tracks[TI.idsml-1].t or \
            Tracks[TI.idsml-1].t==0) or Tracks[TI.idlrg-1].t==0):
        cross,xi,yi =  cross,xi,yi = \
            intersectionTrackTrack(Tracks[TI.idsml-1],Tracks[TI.idlrg-1])
        if cross:
            Tracks[TI.idlrg-1].t = TI.tlrg
            Tracks[TI.idlrg-1].bisect.p2[0] = TI.xi
            Tracks[TI.idlrg-1].bisect.p2[1] = TI.yi
        
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
     
