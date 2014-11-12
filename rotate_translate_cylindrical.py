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

#from __future__ import print_function
import sys
import argparse
import math


def rotate_xaxis(ref_point, orig_point):
  theta = math.atan2(ref_point[2], ref_point[0])
  rotated_point = orig_point
  rotated_point[0] = orig_point[0]*math.cos(math.pi/2 - theta) - orig_point[2]*math.sin(math.pi/2 - theta)
  rotated_point[2] = orig_point[0]*math.sin(math.pi/2 - theta) + orig_point[2]*math.cos(math.pi/2 - theta)
  return rotated_point


def rotate_yaxis(ref_point, orig_point):
  theta = math.atan2(ref_point[2], ref_point[1])
  rotated_point = orig_point
  rotated_point[1] = orig_point[1]*math.cos(math.pi/2 - theta) - orig_point[2]*math.sin(math.pi/2 - theta)
  rotated_point[2] = orig_point[1]*math.sin(math.pi/2 - theta) + orig_point[2]*math.cos(math.pi/2 - theta)
  return rotated_point

def rotate_zaxis(ref_point, orig_point):
  theta = math.atan2(ref_point[1], ref_point[0])
  rotated_point = orig_point
  rotated_point[0] = orig_point[0]*math.cos(-theta) - orig_point[1]*math.sin(-theta)
  rotated_point[1] = orig_point[0]*math.sin(-theta) + orig_point[1]*math.cos(-theta)
  return rotated_point


def xyz_to_cylindrical(xyz):
  cylindrical = []
  cylindrical.append(math.sqrt(pow(float(xyz[0]), 2) + pow(float(xyz[1]), 2)))
  cylindrical.append(math.atan2(float(xyz[1]), float(xyz[0])))
  cylindrical.append(xyz[2])
  return cylindrical


def translate_point(point_to_translate, point_to_subtract):
  translated_point = []
  for x, y in zip(point_to_translate, point_to_subtract):
    translated_point.append(float(x) - float(y))
  return translated_point


def find_max_rad_point(points):
  max_rad_point = points[0]
  max_rad = 0
  for point in points:
    rad = math.sqrt(pow(point[0], 2) + pow(point[1], 2))
    if rad > max_rad:
      max_rad_point = point
      max_rad = rad
  
  return max_rad_point


def read_points(point_file):
  # initialize blank list of lists
  points = []
  for line in point_file :
    point = line.split('\t')
    points.append(point[0:3])
  return points


def write_points(point_file, orig_points, transformed_points):
  assert(len(orig_points) == len(transformed_points))
  point_file.write("orig_x\torig_y\torig_z\trotated_x\trotated_y\trotated_z\trotated_rad\trotated_theta\n")
  for i in range(len(orig_points)):
    orig_point = orig_points[i]
    transformed_point = transformed_points[i]
    cylindrical_point = xyz_to_cylindrical(transformed_point)
    point_file.write("%s\t%s\t%s\n" % ('\t'.join([str(f) for f in orig_point]), '\t'.join([str(f) for f in transformed_point]), '\t'.join([str(f) for f in cylindrical_point[0:2] ])))


def write_point(point_file, orig_point, transformed_point):
  cylindrical_point = xyz_to_cylindrical(transformed_point)
  point_file.write("orig_x\torig_y\torig_z\trotated_x\trotated_y\trotated_z\trotated_rad\trotated_theta\n")
  point_file.write("%s\t%s\t%s\n" % ('\t'.join([str(f) for f in orig_point]), '\t'.join([str(f) for f in transformed_point]), '\t'.join([str(f) for f in cylindrical_point[1:3] ])))
  #point_file.write("%s\t%s\t%s\n" % ('\t'.join(orig_point), '\t'.join(cylindrical_point), '\t'.join(transformed_point)))


def main():
  parser = argparse.ArgumentParser(description='rotate and transform bone scan points')
  parser.add_argument('-p', '--points_infile', dest='points_infile', 
                      type=argparse.FileType('r'),
                      help='Name of input file of all points')
  parser.add_argument('-g', '--origin_infile', dest='origin_infile', 
                      type=argparse.FileType('r'),
                      help='Name of infile corresponding to point to be origin')
  parser.add_argument('-i', '--z_axis_infile', dest='z_axis_infile',
                      type=argparse.FileType('r'),
                      help='Name of infile corresponding to point to be z_axis')
  parser.add_argument('-r', '--rot_ref_infile', dest='rot_ref_infile', 
                      type=argparse.FileType('r'), nargs='?',
                      help='Name of infile corresponding to point to rotate around, if not included the max radius point will be used')
  parser.add_argument('-z', '--z_axis_outfile', dest='z_axis_outfile',
                      type=argparse.FileType('w'),
                      help='Name of outfile corresponding to point to be z_axis')
  parser.add_argument('-o', '--points_outfile', dest='points_outfile', nargs='?',
                      type=argparse.FileType('w'), default=sys.stdout, 
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
  write_point(args.z_axis_outfile, original_zaxis, rotated_zaxis)
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
  #if not rot_ref_infile:
  ref_point = find_max_rad_point(rotated_points)
  #else:
    #temp_points = read_points(args.rot_ref_infile)
    #input_ref_point = temp_points[0]
    #rotated_ref_point = translate_point(input_ref_point, origin)
    #rotated_ref_point = rotate_yaxis(z_axis, rotated_ref_point)
    #ref_point = rotate_xaxis(yrotated_zaxis, rotated_ref_point)
  # check that max_rad_point is set
  assert (len(ref_point) == 3)
  # rotate points so that ref_point is on the y-axis, output points
  full_transformed_points = []
  for rotated_point in rotated_points:
    full_transformed_points.append(rotate_zaxis(ref_point, rotated_point))
  write_points(args.points_outfile, orig_points, full_transformed_points)

if __name__ == '__main__':
  main()
    
    

