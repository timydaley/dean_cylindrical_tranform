# rotate_translate_cylindrical:
#   rotate and translate points from bone scan so that the top
#   and bottom lie on the z-axis and the point furthest away
#   lies on y-axis
#
#    Copyright (C) 2012 University of Southern California and
#			 Timothy Daley
#
#    Authors: Timothy Daley
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.	If not, see <http://www.gnu.org/licenses/>.

import sys
import argparse

def rotate_xaxis(ref_point, orig_point):
  theta = atan2(ref_point[2], ref_point[0])
  rotated_point = orig_point
  rotated_point[0] = orig_point[0]*cos(pi/2 - theta) - orig_point[2]*sin(pi/2 - theta)
  rotated_point[2] = orig_point[0]*sin(pi/2 - theta) + orig_point[2]*cos(pi/2 - theta)
  return rotated_point


def rotate_yaxis(ref_point, orig_point):
  theta = atan2(ref_point[2], ref_point[1])
  rotated_point = orig_point
  rotated_point[1] = orig_point[1]*cos(pi/2 - theta) - orig_point[2]*sin(pi/2 - theta)
  rotated_point[2] = orig_point[1]*sin(pi/2 - theta) + orig_point[2]*cos(pi/2 - theta)
  return rotated_point

def rotate_zaxis(ref_point, orig_point):
  theta = atan2(ref_point[1], ref_point[0])
  rotated_point = orig_point
  rotated_point[0] = orig_point[0]*cos(-theta) - orig_point[1]*sin(-theta)
  rotated_point[1] = orig_point[0]*sin(-theta) + orig_point[1]*cos(-theta)
  return rotated_point


def xyz_to_cylindrical(xyz):
  cylindrical = xyz
  cylindrical[0] = sqrt(pow(xyz[0], 2) + pow(xyz[1], 2))
  cylindrical[1] = atan2(xyz[1], xyz[0])
  return cylindrical


def translate_point(point_to_translate, point_to_subtract):
  translated_point = []
  for x, y in zip(point_to_translate, point_to_subtract):
    translated_point = x - y
  return translated_point


def find_max_rad_point(points):
  max_rad_point = points[0]
  max_rad = 0
  for point in points:
    rad = sqrt(pow(point[0], 2) + pow(point[1], 2))
    if rad > max_rad:
      max_rad_point = point
      max_rad = rad
  return max_rad_point


def read_points(filename):
  # initialize blank list of lists
  points = []
  with open(filename) as file:
    for line in file :
      point = line.split('\t')
      points.append(point)
  return points


def write_points(filename, orig_point, transformed_points):
  file = open(filename, 'w')
  for orig_point, transformed_point in zip(orig_points, transformed_points): 
    cylindrical_point = xyz_to_cylindrical(orig_point)
    file.write("%s\t%s\t%s\n" % (orig_point, cylindrical_point, transformed_point))


def main():
  parser = argparse.ArgumentParser(description='rotate and transform bone scan points')
  parser.add_argument('--points_infile', dest='points_infile', 
                      metavar='p', nargs='+',
                      help='Name of input file of all points')
  parser.add_argument('--origin_infile', dest='origin_infile', 
                      metavar='g', nargs='+',
                      help='Name of infile corresponding to point to be origin')
  parser.add_argument('--z_axis_infile', dest='z_axis_infile', 
                      metavar='i', nargs='+',
                      help='Name of infile corresponding to point to be z_axis')
  parser.add_argument('--rot_ref_infile', dest='rot_ref_infile', 
                      metavar='r', nargs='+',
                      help='Name of infile corresponding to point to rotate around, if not included the max radius point will be used')
  parser.add_argument('--z_axis_outfile', dest='z_axis_outfile', 
                      metavar='z', nargs='+',
                      help='Name of outfile corresponding to point to be z_axis')
  parser.add_argument('--rot_ref_outfile', dest='rot_ref_outfile', 
                      metavar='f', nargs='+',
                      help='Name of outfile for rotated reference point')
  parser.add_argument('--points_outfile', dest='points_outfile', 
                      metavar='o', nargs='+', default=sys.stdout, 
                      help='Name of output file (default: stdout)')
  args = parser.parse_args()
  # read in origin
  temp_points = read_points(args.origin_infile)
  origin = temp_points[0]
  # check length of origin
  assert (len(origin) == 3)
  # read in future z-axis
  temp_points = read_points(args.z_axis_infile)
  z_axis = temp_points[0]
  # check length of z-axis
  assert (len(z_axis) == 3)
  original_zaxis = z_axis
  # translate future z-axis (origin will be (0,0,0) )
  z_axis = translate_point(z_axis, origin)
  # rotate future z-axis to make it present z-axis
  yrotated_zaxis = rotate_yaxis(z_axis, z_axis)
  rotated_zaxis = rotate_xaxis(yrotated_zaxis, yrotated_zaxis)
  # write z-axis
  write_points(args.z_axis_outfile, origonal_zaxis, rotated_zaxis)
  # now read in the rest of the points and transform them
  orig_points = read_points(args.points_infile)
  rotated_points = []
  for point in orig_points:
    # check length of each point
    assert (len(point) == 3)
    # translate point
    rotated_point = point
    rotated_point = translate_point(rotated_point, origin)
    # rotate point add it to rotated_points
    rotated_point = rotate_yaxis(z_axis, rotated_point)
    rotated_points.append(rotate_xaxis(yrotated_zaxis, rotated_point))
   # check if rot_ref_infile was specified, if not find it by max radius
  ref_point = []
  if not rot_ref_infile:
    ref_point = find_max_rad_point(rotated_points)
  else:
    temp_points = read_points(args.rot_ref_infile)
    input_ref_point = temp_points[0]
    rotated_ref_point = translate_point(input_ref_point, origin)
    rotated_ref_point = rotate_yaxis(z_axis, rotated_ref_point)
    ref_point = rotate_xaxis(yrotated_zaxis, rotated_ref_point)
  # check that max_rad_point is set
  assert (len(ref_point) == 3)
  # rotate points so that ref_point is on the y-axis, output points
  full_transformed_points = []
  for rotated_point in rotated_points:
    full_transformed_points.append(rotate_zaxis(ref_point, rotated_point))
  write_points(orig_points, transformed_points)

if __name__ == '__main__':
  main()
    
    

