"""Exact bounded periodic search using the generated radius-1 collar SVG metadata.
No SAT dependency: bitset arc consistency plus exhaustive finite-domain search.
Every solution is checked again directly against the 14 rooted collars.
"""
import argparse
import collections
import json
import pathlib
import time
import xml.etree.ElementTree as ET

DIRS=((1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1))
OFFSETS=((0,0),)+DIRS
COLORS=('yellow','blue','green','purple','red')
NS='{http://www.w3.org/2000/svg}'

def variants(path):
    root=ET.parse(path).getroot()
    assert root.attrib['data-radius']=='1' and root.attrib['data-count']=='14'
    result=set()
    panels=[p for p in root.iter(NS+'g') if p.attrib.get('class')=='collar-panel']
    assert len(panels)==14
    for panel in panels:
        raw={}
        for tile in panel:
            if tile.attrib.get('class')!='hex-tile': continue
            a=tile.attrib
            raw[int(a['data-q']),int(a['data-r'])]=(COLORS.index(a['data-color']),int(a['data-orientation']))
        assert set(raw)==set(OFFSETS)
        for turn in range(6):
            turned={}
            for (q,r),(color,orientation) in raw.items():
                for _ in range(turn): q,r=-r,q+r
                period=2 if color in (0,4) else 6
                turned[q,r]=(color,(orientation+turn)%period)
            result.add(tuple(turned[p] for p in OFFSETS))
    return sorted(result)

def bits(mask):
    while mask:
        bit=mask & -mask
        yield bit.bit_length()-1
        mask-=bit

class Limit(Exception): pass

class Search:
    def __init__(self,states,node_limit=100000,time_limit=5):
        self.states=states
        self.full=(1<<len(states))-1
        self.colors=[sum(1<<i for i,s in enumerate(states) if s[0][0]==c) for c in range(5)]
        self.red0=sum(1<<i for i,s in enumerate(states) if s[0]==(4,0))
        self.compat=[]
        for dq,dr in DIRS:
            shared=[(i,OFFSETS.index((q-dq,r-dr))) for i,(q,r) in enumerate(OFFSETS) if (q-dq,r-dr) in OFFSETS]
            self.compat.append([sum(1<<j for j,t in enumerate(states) if all(s[u]==t[v] for u,v in shared)) for s in states])
        self.node_limit=node_limit;self.time_limit=time_limit

    def solve(self,a,b,skew,all_colors=True):
        self.nodes=0;self.deadline=time.monotonic()+self.time_limit
        def index(q,r):
            k,r=divmod(r,b)
            return r*a+(q-k*skew)%a
        self.neighbors=[tuple(index(q+dq,r+dr) for dq,dr in DIRS) for r in range(b) for q in range(a)]
        domains=[self.full]*(a*b)
        # Any all-color periodic tiling can be translated and properly rotated
        # so a red tile with canonical orientation zero is at the origin.
        if all_colors: domains[0]=self.red0
        self.all_colors=all_colors
        def propagate(ds):
            queue=collections.deque(range(a*b));pending=set(queue)
            while True:
                while queue:
                    x=queue.popleft();pending.remove(x)
                    for direction,y in enumerate(self.neighbors[x]):
                        table=self.compat[direction]
                        if x==y:
                            revised=sum(1<<s for s in bits(ds[x]) if table[s] & (1<<s))
                        else:
                            revised=sum(1<<s for s in bits(ds[x]) if table[s] & ds[y])
                        if not revised: return False
                        if revised!=ds[x]:
                            ds[x]=revised
                            for z in self.neighbors[x]:
                                if z not in pending: queue.append(z);pending.add(z)
                changed=False
                if self.all_colors:
                    for color_mask in self.colors:
                        possible=[i for i,d in enumerate(ds) if d & color_mask]
                        if not possible: return False
                        if len(possible)==1:
                            x=possible[0];revised=ds[x] & color_mask
                            if revised!=ds[x]:
                                ds[x]=revised;changed=True
                                for z in self.neighbors[x]:
                                    if z not in pending:queue.append(z);pending.add(z)
                if not changed:return True
        def visit(ds):
            self.nodes+=1
            if self.nodes>self.node_limit or time.monotonic()>self.deadline:raise Limit
            if not propagate(ds):return None
            undecided=[i for i,d in enumerate(ds) if d.bit_count()>1]
            if not undecided:return [d.bit_length()-1 for d in ds]
            x=min(undecided,key=lambda i:ds[i].bit_count())
            for state in bits(ds[x]):
                child=ds.copy();child[x]=1<<state
                found=visit(child)
                if found is not None:return found
            return None
        try: solution=visit(domains)
        except Limit:return 'unknown',None,self.nodes
        if solution is None:return 'unsat',None,self.nodes
        labels=[self.states[s][0] for s in solution]
        # Independent validation of the wrapped neighborhood at each site,
        # ignoring the compatibility tables used by the search.
        allowed=set(self.states)
        for r in range(b):
            for q in range(a):
                neighborhood=tuple(labels[index(q+dq,r+dr)] for dq,dr in OFFSETS)
                assert neighborhood in allowed
        if all_colors:assert {c for c,o in labels}==set(range(5))
        return 'sat',labels,self.nodes

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--catalog',default='out/build/collars-release/collar-test-output/collars-1.svg')
    p.add_argument('--area',type=int,default=100)
    p.add_argument('--seconds',type=float,default=2)
    p.add_argument('--nodes',type=int,default=100000)
    p.add_argument('--output',default='out/build/radius1-periodic.json')
    args=p.parse_args()
    states=variants(args.catalog);search=Search(states,args.nodes,args.seconds)
    print('oriented radius-1 states',len(states),flush=True)
    checked=0;unknown=[]
    for area in range(5,args.area+1):
        for a in range(1,area+1):
            if area%a:continue
            b=area//a
            for skew in range(a):
                status,labels,nodes=search.solve(a,b,skew)
                checked+=1
                if status=='unknown':unknown.append((a,b,skew))
                if status=='sat':
                    result=dict(radius=1,periods=[[a,0],[skew,b]],area=area,
                                cells=[dict(q=q,r=r,color=COLORS[labels[r*a+q][0]],orientation=labels[r*a+q][1]) for r in range(b) for q in range(a)],
                                checked=checked,unresolved=unknown)
                    pathlib.Path(args.output).write_text(json.dumps(result,indent=2))
                    print('FOUND',result['periods'],'area',area,'counts',dict(collections.Counter(x['color'] for x in result['cells'])),'checked',checked,'unresolved',len(unknown),flush=True)
                    return
        print('area',area,'completed; cells tested',checked,'unresolved',len(unknown),flush=True)
    print('No witness through area',args.area,'unresolved cases',unknown,flush=True)

if __name__=='__main__':main()
