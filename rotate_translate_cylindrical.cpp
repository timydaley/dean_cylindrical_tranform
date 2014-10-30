/*    rotate_translate_cylindrical:
 *
 *    Copyright (C) 2012 University of Southern California and
 *			 Timothy Daley
 *
 *    Authors: Timothy Daley
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.	If not, see <http://www.gnu.org/licenses/>.
 */
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"

#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>
#include <sstream>
#include <math.h>

#define _USE_MATH_DEFINES


using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;

using std::fabs;
using std::istringstream;

static const char row_delim = '\n';
static const char field_delim = '\t';


static void
input_point(const string infile,
	    vector<double> &points){

  points.clear();
  std::ifstream in(infile.c_str());
  if(!in)
    throw SMITHLABException("problem opening input file " 
			   + infile);
  string in_row;
  getline(in, in_row, row_delim);
  istringstream row_stream(in_row);
  for(string field; getline(row_stream, field, field_delim);)
    points.push_back(atof(field.c_str()));
  if(points.size() != 3)
    throw SMITHLABException("problem opening input file " 
			   + infile);
}


//rotation about y-axis
static void
rotate_xaxis(const vector<double> &reference_point,
	     const vector<double> &orig_point,
	     vector<double> &rotated_point){
  //compute angle on x-z plane
  const double theta = atan2(reference_point[2], reference_point[0]);
  rotated_point = orig_point;
  rotated_point[0] = 
    orig_point[0]*cos(M_PI/2.0 - theta) - orig_point[2]*sin(M_PI/2.0 - theta);
  rotated_point[2] = 
    orig_point[0]*sin(M_PI/2.0 - theta) + orig_point[2]*cos(M_PI/2.0 - theta);
}

//rotation about x-axis
static void
rotate_yaxis(const vector<double> &reference_point, 
	     const vector<double> &orig_point,
	     vector<double> &rotated_point){
  //compute angle on y-z plane
  const double theta = atan2(reference_point[2], reference_point[1]);
  rotated_point = orig_point;
  rotated_point[1] = 
    orig_point[1]*cos(M_PI/2.0 - theta) - orig_point[2]*sin(M_PI/2.0 - theta);
  rotated_point[2] = 
    orig_point[1]*sin(M_PI/2.0 - theta) + orig_point[2]*cos(M_PI/2.0 - theta);
}

static void
rotate_zaxis(const vector<double> &reference_point,
	     const vector<double> &orig_point,
	     vector<double> &rotated_point){
  //compute angle on y-z plane
  const double theta = atan2(reference_point[1], reference_point[0]);
  rotated_point.resize(3, 0.0);
  rotated_point[2] = orig_point[2];
  rotated_point[0] = 
    orig_point[0]*cos(-theta) - orig_point[1]*sin(-theta);
  rotated_point[1] = 
    orig_point[0]*sin(-theta) + orig_point[1]*cos(-theta);
}


/*
//rotate about x&y axises
static void
rotate_xyaxis(const vector<double> &zaxis_point,
	      const vector<double> &orig_point,
	      vector<double> &rotated_point){
  const double y_theta = M_PI/2.0 - atan2(zaxis_point[2], zaxis_point[0]);
  const double x_theta = M_PI/2.0 - atan2(zaxis_point[2], zaxis_point[1]);

  rotated_point.resize(3, 0.0);
  rotated_point[0] =
    orig_point[0]*cos(y_theta) + orig_point[1]*sin(x_theta)*sin(y_theta)
    + orig_point[2]*cos(x_theta)*sin(y_theta);
  rotated_point[1] = 
    orig_point[1]*cos(x_theta) - orig_point[2]*sin(x_theta);
  rotated_point[2] = 
    -orig_point[0]*sin(y_theta) + orig_point[1]*sin(x_theta)*cos(y_theta)
    + orig_point[2]*cos(x_theta)*cos(y_theta);
}
*/

/*
// spherical_point[0] = radius
// spherical_point[1] = azimuthal angle(angle of point projected on xy plane)
// spherical_point[2] = polar angle(angle of point projected on xz plane)
static void
xyz_to_spherical(const vector<double> &xyz_point,
		 vector<double> &spherical_point){
  spherical_point.resize(3, 0.0);
  spherical_point[0] = sqrt(xyz_point[0]*xyz_point[0]
			    + xyz_point[1]*xyz_point[1]
			    + xyz_point[2]*xyz_point[2]);
  spherical_point[1] = atan2(xyz_point[1], xyz_point[0]);
  if(spherical_point[0] > 0)
    spherical_point[2] = acos(xyz_point[2]/spherical_point[0]);
}

//cylindrical_point[0] = radius of point(of the point projected on xy plane)
//cylindrical_point[1] = angle of point(projected on xy plane)
//cylindrical_point[2] = z
static void
spherical_to_cylindrical(const vector<double> &spherical_point,
			vector<double> &cylindrical_point){
  cylindrical_point.resize(3, 0.0);
  cylindrical_point[0] = fabs(spherical_point[0]*sin(spherical_point[2]));
  cylindrical_point[1] = spherical_point[1];
  cylindrical_point[2] = spherical_point[0]*cos(spherical_point[2]);
}

static void
cylindrical_to_xyz(const vector<double> &cylindrical_point,
		   vector<double> &xyz_point){
  xyz_point.resize(3, 0.0);
  xyz_point[2] = cylindrical_point[2];
  xyz_point[0] = cylindrical_point[0]*cos(cylindrical_point[1]);
  xyz_point[1] = cylindrical_point[0]*sin(cylindrical_point[1]);
}
*/

static void
xyz_to_cylindrical(const vector<double> &xyz_point,
		   vector<double> &cylindrical_point){
  cylindrical_point.resize(3, 0.0);
  cylindrical_point[2] = xyz_point[2];
  cylindrical_point[0] = 
    sqrt(pow(xyz_point[0], 2) + pow(xyz_point[1], 2));
  cylindrical_point[1] = atan2(xyz_point[1], xyz_point[0]);
}

/*
static void
spherical_to_xyz(const vector<double> &spherical_point,
		 vector<double> &xyz_point){
  xyz_point.resize(3, 0.0);
  xyz_point[0] = 
    spherical_point[0]*sin(spherical_point[2])*cos(spherical_point[1]);
  xyz_point[1] = 
    spherical_point[0]*sin(spherical_point[2])*sin(spherical_point[1]);
  xyz_point[2] = 
    spherical_point[0]*cos(spherical_point[2]);
}
*/

// rotate all points so that point with maximum radius
// is on x-axis
static size_t
find_max_radius(const bool VERBOSE,
		const vector< vector<double> > &points){

  size_t max_rad_indx;
  // find max radius and corresponding point
  double max_rad = 0.0;
  for(size_t i = 0; i < points.size(); i++){
    const double radius = 
      sqrt(pow(points[i][0], 2) + pow(points[i][1], 2));
    if(radius > max_rad){
      max_rad = radius;
      max_rad_indx = i;
    }
  }
  return max_rad_indx;
}



static void
write_point(const string outfile,
	    const vector<double> &orig_xyz_point,
	    const vector<double> &rotated_point,
	    const size_t out_precision){
  std::ofstream of;
  if(!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty()? std::cout.rdbuf() : of.rdbuf());

  out.setf(std::ios_base::fixed, std::ios_base::floatfield);
  out.precision(out_precision);

  vector<double> rotated_cyl_point(3, 0.0);
  xyz_to_cylindrical(rotated_point, rotated_cyl_point);
  out << "original_x" << "\t" << "orginal_y" << "\t" << "original_z" << "\t"
	<< "flipped_x" << "\t" << "flipped_y" << "\t" << "flipped_z" 
	<< "\t" << "flipped_radius" << "\t" << "flipped_theta" <<  endl;
  out << orig_xyz_point[0] << "\t" << orig_xyz_point[1] << "\t" 
      << orig_xyz_point[2] << "\t" << rotated_point[0] << "\t"
      << rotated_point[1] << "\t" << rotated_point[2] << "\t"
      << rotated_cyl_point[0] << "\t" << rotated_cyl_point[1] << endl;
}



int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile, origin_infile, zaxis_infile, ref_infile, origin_outfile, zaxis_outfile, ref_outfile;
    bool VERBOSE = false;
    size_t out_precision = 4;


    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("rotate_translate_cylindrical", "",
			   "<bed-format-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("origin_in", 'g', 
		      "Name of infile corresponding to point to be origin",
		      false, origin_infile);
    opt_parse.add_opt("z_axis_in", 'i', 
		      "Name of infile corresponding to point to be z_axis",
		      false, zaxis_infile);
    opt_parse.add_opt("rot_ref_in",'r',
		      "Name of infile corresponding to point to rotate around, if not included the max radius point will be used",
		      false, ref_infile);
    opt_parse.add_opt("z_axis_out", 'z', 
		      "Name of outfile corresponding to point to be z_axis",
		      false, zaxis_outfile);
    opt_parse.add_opt("rot_ref_out",'f',
		      "Name of outfile for rotated reference point",
		      false, ref_outfile);
    opt_parse.add_opt("out_prec",'p', "number of digits to print in output, default 4", false, out_precision);
    opt_parse.add_opt("verbose", 'v', "print more run information", 
		      false , VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string input_file_name = leftover_args.front();
    /**********************************************************************/
    

    // input origin
    vector<double> origin;
    input_point(origin_infile, origin);
    if(VERBOSE)
      cerr << "origin\t" << origin[0] << "\t" << origin[1] << "\t" << origin[2] << endl;


    // input z_axis
    vector<double> z_axis, orig_zaxis;
    input_point(zaxis_infile, z_axis);
    orig_zaxis = z_axis;
    for(size_t i = 0; i < z_axis.size(); i++)
      z_axis[i] -= origin[i];

    //make suere zaxis gets rotated correctly
    vector<double> rotated_zaxis, yrotated_zaxis;
    rotate_yaxis(z_axis, z_axis, yrotated_zaxis);
    rotate_xaxis(yrotated_zaxis, yrotated_zaxis, rotated_zaxis);
    //rotate_xyaxis(z_axis, z_axis, rotated_zaxis);
    if(VERBOSE){
      cerr << "original z-axis: \t"
	   << orig_zaxis[0] << "\t" << orig_zaxis[1] << "\t" << orig_zaxis[2] << endl;
      cerr << "translated z-axis:\t" 
	   << z_axis[0] << "\t" << z_axis[1] << "\t" << z_axis[2] << endl;
      cerr << "y-rotated z-axis: \t" << yrotated_zaxis[0] << "\t" 
	   << yrotated_zaxis[1] << "\t" << yrotated_zaxis[2] << endl;
      cerr << "rotated z-axis: \t" 
	   << rotated_zaxis[0] << "\t" << rotated_zaxis[1]
	   << "\t" << rotated_zaxis[2] << endl;
    }

    write_point(zaxis_outfile, orig_zaxis, rotated_zaxis, out_precision);

    std::ifstream in(input_file_name.c_str());
    if(!in)
      throw SMITHLABException("problem opening input file " 
			     + input_file_name);
    vector< vector<double> > original_points, translated_points, flipped_points;
    for(string in_row; getline(in, in_row, row_delim); ){
      vector<double> point;
      istringstream row_stream(in_row);
      for(string field; getline(row_stream, field, field_delim);)
	point.push_back(atof(field.c_str()));
      if(point.size() != 3){
	if(point.back() == 0)
	  point.pop_back();
	if(point.size() != 3){
	  cerr << "too few dimensions at line " << original_points.size() + 1 << endl; 
	  cerr << "point = ";
	  for(size_t i = 0; i < point.size(); i++)
	    cerr << point[i] << ", ";
	  cerr << endl;
	}
      }
      assert(point.size() == 3);
      original_points.push_back(point);
      for(size_t j = 0; j < point.size(); j++)
	point[j] -= origin[j];
      translated_points.push_back(point);
      vector<double> rotated_point = point;
      vector<double> yrotated_point = point;
      rotate_yaxis(z_axis, point, yrotated_point);
      rotate_xaxis(yrotated_zaxis, yrotated_point, rotated_point);
      //rotate_xyaxis(z_axis, point, rotated_point);
      flipped_points.push_back(rotated_point);
    }

    vector< vector<double> > rotated_points;
    if(ref_infile.empty()){
      const size_t max_rad_indx = find_max_radius(VERBOSE, flipped_points);
      const vector<double> ref_point = flipped_points[max_rad_indx];
      if(VERBOSE){
	cerr << "max radius = " << sqrt(pow(ref_point[0], 2) + pow(ref_point[1], 2)) << endl;
	cerr << "acheived @ " << max_rad_indx + 1 << "th point" << endl;
	cerr << "original ref point = (" << original_points[max_rad_indx][0] << "\t"
	     << original_points[max_rad_indx][1] << "\t" << original_points[max_rad_indx][2] << ")"
	     << endl;
      }
      if(VERBOSE){
	cerr << "flipped ref point = (" << flipped_points[max_rad_indx][0] 
	     << "\t" << flipped_points[max_rad_indx][1] << "\t" 
	     << flipped_points[max_rad_indx][2] << ")" << endl;
	cerr << "ref point = (" << ref_point[0] << "\t"
	     << ref_point[1] << "\t" << ref_point[2] << ")" << endl;
      }
      for(size_t i = 0; i < flipped_points.size(); i++){
	vector<double> rotated_point;
	rotate_zaxis(ref_point, flipped_points[i], rotated_point);
	rotated_points.push_back(rotated_point);
      }
      vector<double> rotated_ref_point;
      rotate_zaxis(ref_point, ref_point, rotated_ref_point);
      write_point(ref_outfile, original_points[max_rad_indx], rotated_ref_point, out_precision);
    }
    else{
      vector<double> input_ref_point, ref_point;
      input_point(ref_infile, input_ref_point);
      if(VERBOSE){
	cerr << "input ref point = (" << input_ref_point[0] 
	     << "\t" << input_ref_point[1] << "\t" << input_ref_point[2] << ")" << endl;
      }
      // reference point is in original form
      // subtract origin
      for(size_t j = 0; j < input_ref_point.size(); j++)
	input_ref_point[j] -= origin[j];
      // flip point
      vector<double> yrotated_ref_point, zrotated_ref_point;
      rotate_yaxis(z_axis, input_ref_point, yrotated_ref_point);
      rotate_xaxis(yrotated_zaxis, yrotated_ref_point, ref_point);
      if(VERBOSE){
	cerr << "flipped ref point = (" << ref_point[0] << "\t"
	     << ref_point[1] << "\t" << ref_point[2] << ")" << endl;
      }
      for(size_t i = 0; i < flipped_points.size(); i++){
	vector<double> rot_point = flipped_points[i];
	rotate_zaxis(ref_point, flipped_points[i], rot_point);
	rotated_points.push_back(rot_point);
      }
      rotate_zaxis(ref_point, ref_point, zrotated_ref_point); 
      write_point(ref_outfile, input_ref_point, zrotated_ref_point, out_precision);
    }

    // output
    std::ofstream of;
    if(!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty()? std::cout.rdbuf() : of.rdbuf());

    out.setf(std::ios_base::fixed, std::ios_base::floatfield);
    out.precision(out_precision);

    out << "orig_x" << "\t" << "orig_y" << "\t" << "orig_z" << "\t"
	<< "flipped_x" << "\t" << "flipped_y" << "\t" << "flipped_z" << "\t"
	<< "flipped_rad" << "\t" << "flipped_theta" << "\t"
	<< "rotated_x" << "\t" << "rotated_y" << "\t" << "rotated_z" << "\t"
	<< "rotated_rad" << "\t" << "rotated_theta" << endl;
    for(size_t i = 0; i < original_points.size(); i++){
      for(size_t j = 0; j < original_points[i].size(); j++)
	out << original_points[i][j] << "\t";
      for(size_t j = 0; j < flipped_points[i].size(); j++)
	out << flipped_points[i][j] << "\t";
      vector<double> flipped_cyl_point(3, 0.0);
      xyz_to_cylindrical(flipped_points[i], flipped_cyl_point);
      out << flipped_cyl_point[0] << "\t" << flipped_cyl_point[1] << "\t";
      for(size_t j = 0; j < rotated_points[i].size(); j++)
	out << rotated_points[i][j] << "\t";
      vector<double> rotated_cylindrical_point(3, 0.0);
      xyz_to_cylindrical(rotated_points[i], rotated_cylindrical_point);
      out << rotated_cylindrical_point[0] << "\t" << rotated_cylindrical_point[1];
      out << endl;
    }


  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
