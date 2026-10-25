"""Optional standard-library XML validation; run after cli_validation.cmake."""
import math
import pathlib
import re
import sys
import xml.etree.ElementTree as ET

ROOT = pathlib.Path(sys.argv[1])
NS = '{http://www.w3.org/2000/svg}'
COLORS = {'yellow': '#ffff00', 'blue': '#3fbcec', 'green': '#3da347',
          'purple': '#6327ac', 'red': '#ff0000'}
FILES = list(ROOT.glob('*/*.svg'))
assert FILES, 'No SVG artifacts'
for path in FILES:
    root = ET.parse(path).getroot()
    assert root.tag == NS + 'svg', path
    width, height = float(root.attrib['width']), float(root.attrib['height'])
    assert 0 < width <= 1341 and 0 < height <= 1341, (path, width, height)
    assert list(map(float, root.attrib['viewBox'].split())) == [0, 0, width, height]
    leaves = [e for e in root.iter(NS + 'g') if e.attrib.get('class') == 'subtile']
    assert leaves, path
    for leaf in leaves:
        polygon = leaf.find(NS + 'polygon')
        assert polygon is not None and polygon.attrib['fill'] == COLORS[leaf.attrib['data-type']], path
        coords = [float(x) for x in re.split(r'[ ,]+', polygon.attrib['points'].strip())]
        assert len(coords) >= 6 and len(coords) % 2 == 0, path
        assert all(math.isfinite(x) for x in coords), path
        assert all(0 <= x <= width for x in coords[::2]), path
        assert all(0 <= y <= height for y in coords[1::2]), path
        # Nonzero signed area; a missing/wrong boundary loop cannot pass this.
        points = list(zip(coords[::2], coords[1::2]))
        area = sum(a[0]*b[1]-a[1]*b[0] for a,b in zip(points,points[1:]+points[:1]))
        assert abs(area) > 0.001, path
        for line in leaf.findall(NS + 'line'):
            assert all(math.isfinite(float(line.attrib[k])) for k in ('x1','y1','x2','y2')), path
    if path.name == 'rules.svg':
        assert len(leaves) == 42, (path, len(leaves))  # five parents + 37 children
        assert {e.attrib['data-type'] for e in leaves} == set(COLORS), path
    if path.parent.name.endswith('-0') and path.name == 'patch.svg':
        assert len(leaves) == 1, path
    if path.parent.name.endswith('-1') and path.name == 'patch.svg' and not path.parent.name.startswith('reference-'):
        seed = path.parent.name.rsplit('-',1)[0]
        assert len(leaves) == dict(yellow=5, blue=11, green=11, purple=8, red=2)[seed], path

def pairs(text):
    values = [float(x) for x in re.split(r'[ ,]+', text.strip())]
    return list(zip(values[::2], values[1::2]))

def distance(a, b):
    return math.hypot(a[0]-b[0], a[1]-b[1])

def on_segment(p, a, b):
    return abs(distance(a,p)+distance(p,b)-distance(a,b)) < 1e-8

def groups(root):
    return [e for e in root.iter(NS+'g') if e.attrib.get('class') == 'subtile']

# Compare actual SVG geometry, not just the affine frames. This catches bad cell
# boundary tracing, lost emphasis, layout errors and reflected arrow offsets.
for seed in COLORS:
    for depth in (0,1,2,4):
        for filename in ('patch.svg', 'rules.svg'):
            actual = ET.parse(ROOT/f'{seed}-{depth}'/filename).getroot()
            reference = ET.parse(ROOT/f'reference-{seed}-{depth}'/filename).getroot()
            for attr in ('width','height'):
                assert abs(float(actual.attrib[attr])-float(reference.attrib[attr])) < 1e-8
            aa, rr = groups(actual), groups(reference)
            assert len(aa) == len(rr)
            for a,r in zip(aa,rr):
                assert a.attrib['data-type'] == r.attrib['data-type']
                av = pairs(a.find(NS+'polygon').attrib['points'])
                rv = pairs(r.find(NS+'polygon').attrib['points'])
                reference_segments = list(zip(rv,rv[1:]+rv[:1]))
                assert all(any(on_segment(p,x,y) for x,y in reference_segments) for p in av)
                assert all(any(distance(p,q)<1e-8 for q in av) for p in rv)
                # Grid paths split long straight strokes into unit edges.
                segments = a.find(NS+'path').attrib['d'].strip().split('M ')[1:]
                actual_segments = [pairs(s.replace('L ','')) for s in segments]
                rv = pairs(r.find(NS+'polyline').attrib['points'])
                reference_segments = list(zip(rv,rv[1:]))
                assert all(any(on_segment(x,u,v) and on_segment(y,u,v)
                               for u,v in reference_segments) for x,y in actual_segments)
                assert abs(sum(distance(x,y) for x,y in actual_segments)-
                           sum(distance(x,y) for x,y in reference_segments)) < 1e-8
                al, rl = a.find(NS+'line'), r.find(NS+'line')
                assert al is not None and rl is not None
                for attr in ('x1','y1','x2','y2'):
                    assert abs(float(al.attrib[attr])-float(rl.attrib[attr])) < 1e-8
print(f'Parsed {len(FILES)} SVG files; outlines, emphasized edges and arrows match geometric references')
