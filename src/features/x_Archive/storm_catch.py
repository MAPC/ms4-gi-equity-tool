from copy import deepcopy

import geopandas as gpd
import networkx as nx
from shapely.geometry import LineString, MultiLineString, Point

from stormcatchments.network import Network

#this is from storm catchments but for some reason
#my system can't find it...
#https://github.com/t-ott/stormcatchments/blob/main/stormcatchments/topology.py#L48


def find_floating_points(net: Network) -> gpd.GeoDataFrame:
  '''
  Find and return any points in a Network that are not snapped to a line vertex. These
  floating points cannot be integrated into networking functionality unless they are
  snapped to a line vertex.

  Parameters
  ----------
  net : Network
    A stormcatchments Network object whose point data will be inspected for floating
    points

  Returns
  -------
  floating_pts : gpd.GeoDataFrame
    A GeoDataFrame of any floating points in net.pts
  '''
  floating_pts = []
  for pt in net.pts.itertuples():
    # Get all segments that pt touches, could be on a vertex or between verticies
    touch_segs = net.segments[net.segments.geometry.touches(pt.geometry)]

    if len(touch_segs) == 0:
      floating_pts.append(pt)
      continue

    # Collect segment coordinates as (x, y) tuples
    seg_coords = set()
    for l in touch_segs.geometry:
      for c in l.coords:
        seg_coords.add(c)

    if (pt.geometry.x, pt.geometry.y) not in seg_coords:
      floating_pts.append(pt)
  
  return gpd.GeoDataFrame(floating_pts, crs=net.crs).set_index('Index')

def snap_points(net: Network, tolerance: float) -> Network:
  '''
  Create a copy of a supplied Network which snaps any points in Network to the
  closest line vertex within a snapping tolerance

  Parameters
  ----------
  net : Network
    A stormcatchments Network object which may contain floating points

  tolerance : float
    The maximum search distance to find the nearest vertex

  Returns
  -------
  net_snapped : Network
    A stormcatchments Network object with snapping applied to its point data
  '''
  # Snap all floating (un-snapped) points to the nearest line vertex
  floating_pts = find_floating_points(net)

  net_snapped = deepcopy(net)
  for pt in floating_pts.itertuples():
      nearby = net.segments.cx[
        pt.geometry.x-tolerance:pt.geometry.x+tolerance,
        pt.geometry.y-tolerance:pt.geometry.y+tolerance
      ]

      closest_xy = None
      closest_dist = tolerance**2
      for l in nearby.geometry:
        for c in l.coords:
          dist = pt.geometry.distance(Point(c))
          if dist < closest_dist:
            closest_dist = dist
            closest_xy = c
      
      if closest_dist <= tolerance:
        net_snapped.pts.at[pt.Index, 'geometry'] = Point(closest_xy)

  return net_snapped