"""Validate the standalone contact sheets as XML and inspect every visible primitive."""
import math
import pathlib
import sys
import xml.etree.ElementTree as ET
NS = '{http://www.w3.org/2000/svg}'
colors = {'yellow': '#ffff00', 'blue': '#3fbcec', 'green': '#3da347', 'purple': '#6327ac', 'red': '#ff0000'}
paths = list(pathlib.Path(sys.argv[1]).glob('*.svg'))
assert {'collars-0.svg', 'collars-1.svg', 'collars-6.svg', 'empty.svg', 'single-column.svg'} <= {p.name for p in paths}
for path in paths:
    root = ET.parse(path).getroot()
    assert root.tag == NS+'svg'
    radius, count = int(root.attrib['data-radius']), int(root.attrib['data-count'])
    expected = 1+3*radius*(radius+1)
    fixed_counts = {'collars-0.svg': 5, 'collars-1.svg': 14, 'collars-6.svg': 149, 'empty.svg': 0, 'single-column.svg': 14}
    if path.name in fixed_counts:
        assert count == fixed_counts[path.name]
    width, height = float(root.attrib['width']), float(root.attrib['height'])
    assert list(map(float,root.attrib['viewBox'].split())) == [0,0,width,height]
    ids=[e.attrib['id'] for e in root.iter() if 'id' in e.attrib]
    assert len(ids)==len(set(ids))
    panels=[e for e in root.iter(NS+'g') if e.attrib.get('class')=='collar-panel']
    assert len(panels)==count
    assert root.find(NS+'title') is not None
    boxes=[]
    for index,panel in enumerate(panels):
        assert int(panel.attrib['data-collar-index'])==index
        label=panel.find(NS+'text')
        assert label.text=='#'+str(index)
        lx,ly=float(label.attrib['x']),float(label.attrib['y'])
        assert 0<lx<width-30 and 14<ly<height
        tiles=[e for e in panel if e.attrib.get('class')=='hex-tile']
        assert len(tiles)==expected
        assert sum(e.attrib['data-central']=='true' for e in tiles)==1
        seen=set();xs=[];ys=[]
        for tile in tiles:
            a=tile.attrib
            q,r=int(a['data-q']),int(a['data-r'])
            assert (q,r) not in seen and max(abs(q),abs(r),abs(q+r))<=radius
            seen.add((q,r))
            assert int(a['data-collar-index'])==index
            orientation=int(a['data-orientation'])
            assert 0<=orientation<(2 if a['data-color'] in ('yellow','red') else 6)
            poly=tile.find(NS+'polygon')
            assert poly.attrib['fill']==colors[a['data-color']]
            points=[tuple(map(float,p.split(','))) for p in poly.attrib['points'].split()]
            assert len(points)==6
            stroke=float(poly.attrib['stroke-width'])/2
            for x,y in points:
                assert math.isfinite(x+y) and stroke<x<width-stroke and stroke<y<height-stroke
                xs.append(x);ys.append(y)
            assert (a['data-central']=='true')==((q,r)==(0,0))
            assert float(poly.attrib['stroke-width'])==(3 if (q,r)==(0,0) else .7)
            line=tile.find(NS+'line')
            assert line.attrib['class']=='orientation-arrow' and line.attrib['marker-end']=='url(#collar-arrowhead)'
            x1,y1,x2,y2=[float(line.attrib[k]) for k in ('x1','y1','x2','y2')]
            cx=sum(x for x,y in points)/6;cy=sum(y for x,y in points)/6
            assert abs(x1-cx)<1e-8 and abs(y1-cy)<1e-8
            angle=math.atan2(-(y2-y1),x2-x1)
            assert abs(math.cos(angle)-math.cos(orientation*math.pi/3))<1e-8
            assert abs(math.sin(angle)-math.sin(orientation*math.pi/3))<1e-8
            # Arrow length .54R plus its .20R marker fits inside the polygon.
            size=math.hypot(points[0][0]-cx,points[0][1]-cy)
            assert math.hypot(x2-cx,y2-cy)+.20*size < math.sqrt(3)/2*size
        assert ly+4<min(ys), 'label overlaps hexes'
        box=(min(xs)-2,min(ys)-2,max(xs)+2,max(ys)+2)
        for x0,y0,x3,y3 in boxes:
            assert box[2]<x0 or box[0]>x3 or box[3]<y0 or box[1]>y3, 'overlapping panels'
        boxes.append(box)
    print(path.name, count, 'panels;', count*expected, 'hexes; XML and bounds valid')

