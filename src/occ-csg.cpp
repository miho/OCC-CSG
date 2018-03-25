/*
 * Copyright 2018 Michael Hoffer <info@michaelhoffer.de>. All rights reserved.
 * 
 * This file is part of OCC-CSG.
 * 
 * OCC-CSG is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 as published by the
 * Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

// std
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string> 
#include <vector>
// #include <filesystem> currently still unusable

// numbers, limits and errors
#include <errno.h>
#include <limits>
#include <stdlib.h>

// primitive objects
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakePrism.hxx>

// CSG operators
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Common.hxx>

// sewing & solid
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <TopoDS_Shell.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>

// STEP import and export
#include <BRepGProp.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <GProp_GProps.hxx>

// BREP import and export
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>

// STL export
#include <StlAPI_Writer.hxx>
#include <BRepMesh_IncrementalMesh.hxx>

// STL import
#include <RWStl.hxx>
#include <StlMesh_Mesh.hxx>
#include <StlMesh_MeshExplorer.hxx>

// relevant for importing stl files
#include <OSD_Path.hxx>

// geometric objects
#include <TopoDS_Compound.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>

#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>

#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>

#include <gp_Circ.hxx>

// fonts & text
#include <Font_BRepFont.hxx>

// TODO we have to upgrade to occt-7.2.x before using this class :(
//#include <Font_BRepTextBuilder.hxx>

// transform
#include <BRepBuilderAPI_GTransform.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_GTrsf.hxx>

// math (OCCT/OCE compliant)
#include <math.hxx>
#define MAX2(X, Y)	(  Abs(X) > Abs(Y)? Abs(X) : Abs(Y) )
#define MAX3(X, Y, Z)	( MAX2 ( MAX2(X,Y) , Z) )

// version
#define VERSION 0.6

// minimal API for primitive objects
TopoDS_Shape createBox(double x1, double y1, double z1, double x2, double y2, double z2);
TopoDS_Shape createSphere(double x1, double y1, double z1, double r);
TopoDS_Shape createCylinder(double r, double h);
TopoDS_Shape createCylinder(double r, double h, double angle);
TopoDS_Shape createCone(double r1, double r2, double h);
TopoDS_Shape createCone(double r1, double r2, double h, double angle);
TopoDS_Shape extrudePolygon(double ex, double ey, double ez, std::vector<double> const &points);
TopoDS_Shape extrudeFile(double ex, double ey, double ez, std::string const &filename);
TopoDS_Shape createCircle(double x, double y, double z, double dx, double dy, double dz, double r);
TopoDS_Shape createPolygon2d(std::vector<double>const &coords);
TopoDS_Shape createRect2d(double minX, double minY, double maxX, double maxY);
TopoDS_Shape createText2d(std::string const &font, double fSize, double x, double y, std::string const& text);

// minimal transform API
TopoDS_Shape transform(TopoDS_Shape shape, double transform_matrix[12]);

// CLI functions
void version();
void usage();
void error();
void notImplemented();
void create(int argc, char *argv[]);
void transform(int argc, char *argv[]);
void convert(int argc, char *argv[]);
void csg(int argc, char *argv[]);
void bounds(int argc, char *argv[]);

// minimal IO API
TopoDS_Shape load(std::string const &filename);
bool save(std::string const &filename, TopoDS_Shape shape, double stlTOL);
TopoDS_Shape importSTL(std::string const &file );
void unsupportedFormat(std::string const &filename);
bool isAccessible(std::string const &filename);

// String API

// checks whether str ends with ending
bool endsWith(std::string const & str, std::string const &ending);

// returns lowercase version of str
std::string toLower(std::string const &str);

// split string by specified separator
std::vector<std::string> split(std::string const &str, const char sep);

// error codes for number conversion
enum NUMBER_CONVERSION_ERROR {
	VALID                      =  0,
	ERR_OUT_OF_RANGE           = -1,
	ERR_CANNOT_PARSE_NUMBER    = -2,
	ERR_EXTRA_CHARS_AT_THE_END = -3
};

// Java-style number conversion (with error checks)
// (atof,ato,... are very buggy and mostly useless, remember to use strto*)
double parseDouble(std::string const &str, NUMBER_CONVERSION_ERROR *ERROR);
double parseDouble(std::string const &str);
int parseInt(std::string const &str, NUMBER_CONVERSION_ERROR *ERROR);
int parseInt(std::string const &str);

// the CLI appliction
int main(int argc, char *argv[])
{
	
	if(argc > 1 && strcmp(argv[1], "--version")==0) { version(); exit(0); }
	
	std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "------        CSG based on the OCE CAD Kernel         ------" << std::endl;
	std::cout << "------                  Version " << VERSION << "                   ------" << std::endl;
	std::cout << "------ 2018 by Michael Hoffer (info@michaelhoffer.de) ------" << std::endl;
	std::cout << "------                www.mihosoft.eu                 ------" << std::endl;
	std::cout << "------------------------------------------------------------" << std::endl;

	if(argc < 2) {
		error();
	}

	if(strcmp(argv[1], "--create")==0) create(argc,argv);
	else if(strcmp(argv[1], "--transform")==0) transform(argc,argv);
	else if(strcmp(argv[1], "--convert")==0) convert(argc,argv);
	else if(strcmp(argv[1], "--csg")==0) csg(argc,argv);
	else if(strcmp(argv[1], "--bounds")==0) bounds(argc,argv);
	else if(strcmp(argv[1], "--help")==0 || strcmp(argv[1], "-h")==0) usage();
	else error();

	return 0;
}

TopoDS_Shape importSTL( std::string const &file )
{
	std::cout << "> importing STL file" << std::endl;

    TCollection_AsciiString aName( (Standard_CString)file.data() );
    OSD_Path aFile(aName);

	BRepBuilderAPI_Sewing shapeSewer;

	std::cout << " -> reading '" << file << "'" << std::endl;

    Handle(StlMesh_Mesh) aSTLMesh = RWStl::ReadFile(aFile);

    Standard_Integer NumberDomains = aSTLMesh->NbDomains();
    
    gp_XYZ p1, p2, p3;
    TopoDS_Vertex Vertex1, Vertex2, Vertex3;
	TopoDS_Shape shape;
    TopoDS_Face face;
    TopoDS_Wire wire;
    Standard_Real x1, y1, z1;
    Standard_Real x2, y2, z2;
    Standard_Real x3, y3, z3;

    StlMesh_MeshExplorer aMExp (aSTLMesh);

	std::cout << " -> converting to faces" << std::endl;

    for (Standard_Integer iND=1;iND<=NumberDomains;iND++)
    {
	for (aMExp.InitTriangle (iND); aMExp.MoreTriangle (); aMExp.NextTriangle ())
	{
	    aMExp.TriangleVertices (x1,y1,z1,x2,y2,z2,x3,y3,z3);
	    p1.SetCoord(x1,y1,z1);
	    p2.SetCoord(x2,y2,z2);
	    p3.SetCoord(x3,y3,z3);

	    if (!p1.IsEqual(p2,0.0) && !p1.IsEqual(p3,0.0))
	    {
			Vertex1 = BRepBuilderAPI_MakeVertex(p1);
			Vertex2 = BRepBuilderAPI_MakeVertex(p2);
			Vertex3 = BRepBuilderAPI_MakeVertex(p3);

			wire = BRepBuilderAPI_MakePolygon( Vertex1, Vertex2, Vertex3, Standard_True);
				if( !wire.IsNull())
				{
					face = BRepBuilderAPI_MakeFace( wire );
					if(!face.IsNull()) {
						shapeSewer.Add(face);
					}
				}
			}
		}
    }

	std::cout << " -> sewing faces" << std::endl;

	shapeSewer.Perform();
	shape = shapeSewer.SewedShape();

	std::cout << " -> extracting shells" << std::endl;
    
	BRepBuilderAPI_MakeSolid solidmaker;
	TopTools_IndexedMapOfShape shellMap;
	TopExp::MapShapes(shape, TopAbs_SHELL, shellMap);

	unsigned int counter = 0;
	for(int ishell = 1; ishell <= shellMap.Extent(); ++ishell) {
    	const TopoDS_Shell& shell = TopoDS::Shell(shellMap(ishell));
    	solidmaker.Add(shell);
		counter++;
	}

	std::cout << "   -> shells found: " << counter << std::endl;

	std::cout << " -> converting to solid" << std::endl;
	
	TopoDS_Shape solid = solidmaker.Solid();

	std::cout << " -> done." << std::endl;

	return solid;
}

TopoDS_Shape load(std::string const &filename) {
	std::cout << "> loading geometry" << std::endl;
	std::cout << " -> reading file '" << filename << "'" << std::endl;

	if(!isAccessible(filename)) {
		std::cerr << "ERROR: file '" << filename << "' cannot be accessed!" << std::endl;
		exit(1);
	}

	TopoDS_Shape shape;

	if(endsWith(toLower(filename), ".stl")) {
		shape = importSTL(filename);
	} else if(endsWith(toLower(filename), ".stp") || endsWith(toLower(filename), ".step")) {
		STEPControl_Reader Reader;
    	Reader.ReadFile(filename.c_str());
    	Reader.TransferRoots();
    	shape = Reader.OneShape();
	} else if(endsWith(toLower(filename), ".brep")){
		BRep_Builder b;
		BRepTools::Read(shape, filename.c_str(), b);
	} else {
		unsupportedFormat(filename);
	}

	std::cout << " -> done." << std::endl;

	return shape;
}

bool save(std::string const &filename, TopoDS_Shape shape, double stlTOL) {
	std::cout << "> saving geometry" << std::endl;
	std::cout << " -> writing file '" << filename << "'" << std::endl;

	if(endsWith(toLower(filename), ".stl")) {
		std::cout << " -> STL TOL: " << stlTOL << std::endl;
		StlAPI_Writer myStlWriter;
		BRepMesh_IncrementalMesh twoMesh( shape, stlTOL);
    	twoMesh.Perform();
		myStlWriter.Write(shape, filename.c_str());
	} else if(endsWith(toLower(filename), ".stp") || endsWith(toLower(filename), ".step")) {
		std::cout << " -> ignoring STL TOL (using resolution independent format): " << stlTOL << std::endl;
		STEPControl_Writer writer;
		writer.Transfer(shape,STEPControl_AsIs);
		writer.Write(filename.c_str());
	}  else if(endsWith(toLower(filename), ".brep")){
		std::cout << " -> ignoring STL TOL (using resolution independent format): " << stlTOL << std::endl;
		BRepTools::Write(shape, filename.c_str());
	} else {
		unsupportedFormat(filename);
	}
	
	std::cout << " -> done." << std::endl;

	return true;
}

void create(int argc, char *argv[]) {

	if(strcmp(argv[2],"box")==0) {
		if(argc != 5 && argc !=6) {
			error();
		}

		std::vector<std::string> values = split(argv[3], ',');

		if(values.size()!=6) {
			error();
		}

		double x1 = parseDouble(values[0].c_str());
		double y1 = parseDouble(values[1].c_str());
		double z1 = parseDouble(values[2].c_str());
		double x2 = parseDouble(values[3].c_str());
		double y2 = parseDouble(values[4].c_str());
		double z2 = parseDouble(values[5].c_str());

		TopoDS_Shape shape = createBox(x1,y1,z1,x2,y2,z2);

		std::string filename = argv[4];

		double stlTOL;

		if(argc == 5) {
			stlTOL = 0.5;
		} else {
			stlTOL = parseDouble(argv[5]);
		}

		save(filename,shape, stlTOL);
	} else if(strcmp(argv[2],"sphere")==0) {
		if(argc != 5 && argc !=6) {
			error();
		}

		std::vector<std::string> values = split(argv[3], ',');

		if(values.size()!=4) {
			error();
		}

		double x1 = parseDouble(values[0].c_str());
		double y1 = parseDouble(values[1].c_str());
		double z1 = parseDouble(values[2].c_str());
		double r  = parseDouble(values[3].c_str());

		TopoDS_Shape shape = createSphere(x1,y1,z1,r);

		std::string filename = argv[4];

		double stlTOL;

		if(argc == 5) {
			stlTOL = 0.5;
		} else {
			stlTOL = parseDouble(argv[5]);
		}

		save(filename,shape, stlTOL);
	} else if(strcmp(argv[2],"cyl")==0) {
		if(argc != 5 && argc !=6) {
			error();
		}

		std::vector<std::string> values = split(argv[3], ',');

		if(values.size()!=5) {
			error();
		}

		double x1 = parseDouble(values[0].c_str());
		double y1 = parseDouble(values[1].c_str());
		double z1 = parseDouble(values[2].c_str());
		double r  = parseDouble(values[3].c_str());
		double h  = parseDouble(values[4].c_str());

		TopoDS_Shape shape = createCylinder(r,h);

		gp_Trsf transformation;
		transformation.SetTranslation(gp_Vec(x1, y1, z1));
		shape = BRepBuilderAPI_GTransform(shape, transformation);

		std::string filename = argv[4];

		double stlTOL;

		if(argc == 5) {
			stlTOL = 0.5;
		} else {
			stlTOL = parseDouble(argv[5]);
		}

		save(filename,shape, stlTOL);
	} else if(strcmp(argv[2],"cone")==0) {
		if(argc != 5 && argc !=6) {
			error();
		}

		std::vector<std::string> values = split(argv[3], ',');

		if(values.size()!=6) {
			error();
		}

		double x1 = parseDouble(values[0].c_str());
		double y1 = parseDouble(values[1].c_str());
		double z1 = parseDouble(values[2].c_str());
		double r1  = parseDouble(values[3].c_str());
		double r2  = parseDouble(values[4].c_str());
		double h  = parseDouble(values[5].c_str());

		TopoDS_Shape shape = createCone(r1,r2,h);

		gp_Trsf transformation;
		transformation.SetTranslation(gp_Vec(x1, y1, z1));
		shape = BRepBuilderAPI_GTransform(shape, transformation);

		std::string filename = argv[4];

		double stlTOL;

		if(argc == 5) {
			stlTOL = 0.5;
		} else {
			stlTOL = parseDouble(argv[5]);
		}

		save(filename,shape, stlTOL);
	} else if(strcmp(argv[2],"extrusion:polygon")==0) {

		if(argc != 5 && argc !=6) {
			error();
		}

		std::vector<std::string> values = split(argv[3], ',');

		std::vector<double> coords;

		// first three entries define the extrusion vector (direction)
		double ex = parseDouble(values[0]);
		double ey = parseDouble(values[1]);
		double ez = parseDouble(values[2]);

		for(size_t i = 3; i < values.size(); i++) {
			double v = parseDouble(values[i]);
			coords.push_back(v);
		}

		TopoDS_Shape shape = extrudePolygon(ex,ey,ez, coords);

		std::string filename = argv[4];

		double stlTOL;

		if(argc == 5) {
			stlTOL = 0.5;
		} else {
			stlTOL = parseDouble(argv[5]);
		}

		save(filename,shape, stlTOL);

	} else if(strcmp(argv[2],"extrusion:file")==0) {

		if(argc != 6 && argc !=7) {
			error();
		}

		std::vector<std::string> values = split(argv[3], ',');
		
		if(values.size()!=3) {
			error();
		}

		// first three entries define the extrusion vector (direction)
		double ex = parseDouble(values[0]);
		double ey = parseDouble(values[1]);
		double ez = parseDouble(values[2]);

		TopoDS_Shape shape = extrudeFile(ex,ey,ez,argv[4]);

		std::string filename = argv[5];

		double stlTOL;

		if(argc == 7) {
			stlTOL = parseDouble(argv[6]);
		} else {
			stlTOL = 0.5;
		}

		save(filename,shape, stlTOL);

	} else if(strcmp(argv[2],"2d:circle")==0) {

		if(argc != 5 && argc !=6) {
			error();
		}

		std::vector<std::string> values = split(argv[3], ',');

		std::vector<double> coords;

		// first three entries define the extrusion vector (direction)
		double x = parseDouble(values[0]);
		double y = parseDouble(values[1]);
		double r = parseDouble(values[2]);

		if(values.size()!=3) {
			error();
		}

		TopoDS_Shape shape = createCircle(x,y,0,0,0,1,r);

		std::string filename = argv[4];

		double stlTOL;

		if(argc == 5) {
			stlTOL = 0.5;
		} else {
			stlTOL = parseDouble(argv[5]);
		}

		save(filename,shape, stlTOL);

	} else if(strcmp(argv[2],"2d:polygon")==0) {

		if(argc != 5 && argc !=6) {
			error();
		}

		std::vector<std::string> values = split(argv[3], ',');

		std::vector<double> coords;

		for(size_t i = 0; i < values.size(); i++) {
			double v = parseDouble(values[i]);
			coords.push_back(v);
		}

		TopoDS_Shape shape = createPolygon2d(coords);

		std::string filename = argv[4];

		double stlTOL;

		if(argc == 5) {
			stlTOL = 0.5;
		} else {
			stlTOL = parseDouble(argv[5]);
		}

		save(filename,shape, stlTOL);

	} else if(strcmp(argv[2],"2d:rect")==0) {

		if(argc != 5 && argc !=6) {
			error();
		}

		std::vector<std::string> values = split(argv[3], ',');

		if(values.size()!=4) {
			error();
		}

		double x1 = parseDouble(values[0]);
		double y1 = parseDouble(values[1]);
		double x2 = parseDouble(values[2]);
		double y2 = parseDouble(values[3]);

		TopoDS_Shape shape = createRect2d(x1,y1,x2,y2);

		std::string filename = argv[4];

		double stlTOL;

		if(argc == 5) {
			stlTOL = 0.5;
		} else {
			stlTOL = parseDouble(argv[5]);
		}

		save(filename,shape, stlTOL);

	} else if(strcmp(argv[2],"2d:text")==0) {

		if(argc != 8 && argc !=9) {
			error();
		}

		std::string fontFileName = argv[3];

		double fontSize = parseDouble(argv[4]);

		std::vector<std::string> values = split(argv[5], ',');

		if(values.size()!=2) {
			error();
		}

		double x = parseDouble(values[0]);
		double y = parseDouble(values[1]);

		std::string text = argv[6];

		TopoDS_Shape shape = createText2d(fontFileName, fontSize, x, y, text);

		std::string filename = argv[7];

		double stlTOL;

		if(argc == 8) {
			stlTOL = 0.5;
		} else {
			stlTOL = parseDouble(argv[8]);
		}

		save(filename,shape, stlTOL);

	} else error();

}

void convert(int argc, char *argv[]) {

  if(argc == 4) {
	  TopoDS_Shape srcShape    = load(argv[2]);
	  save(argv[3], srcShape, 0.5);
  } else if(argc == 5) {
	  TopoDS_Shape srcShape    = load(argv[2]);
	  save(argv[3], srcShape, parseDouble(argv[4]));
  } else {
	  error();
  }

}

void csg(int argc, char *argv[]) {
	if(argc != 6 && argc != 7) {
		error();
	}

	TopoDS_Shape res;

	TopoDS_Shape s1 = load(argv[3]);
	TopoDS_Shape s2 = load(argv[4]);

	std::cout << "> applying csg operation" << std::endl;

	if(strcmp(argv[2],"union")==0) {
		BRepAlgoAPI_Fuse csg(s1, s2);
		res = csg.Shape();
	} else if(strcmp(argv[2],"difference")==0) {
		BRepAlgoAPI_Cut csg(s1, s2);
		res = csg.Shape();
	} else if(strcmp(argv[2],"intersection")==0) {
		BRepAlgoAPI_Common csg(s1, s2);
		res = csg.Shape();
	} else {
		error();
	}

	std::cout << " -> done." << std::endl;

	double stlTOL;

	if(argc == 7) {
		stlTOL = parseDouble(argv[6]);
	} else {
		stlTOL = 0.5;
	}

	save(argv[5], res, stlTOL);
}

void bounds(int argc, char *argv[]) {
	if(argc < 3 || argc > 5) {
		error();
	}

	TopoDS_Shape shape = load(argv[2]);

	std::cout << "> computing bounding box" << std::endl;

	std::cout << " -> approximating bounds" << std::endl;

    // compute bbox on geometric object
	double xMin,yMin,zMin,xMax,yMax, zMax = 0;
	Standard_Real aDeflection = 0.001, deflection;
	Bnd_Box box;
	BRepBndLib::Add(shape, box);
	box.Get(xMin, yMin, zMin, xMax, yMax, zMax);
	deflection= MAX3( xMax-xMin , yMax-yMin , zMax-zMin)*aDeflection;
	
	std::cout << " -> tesselating object" << std::endl;
	// tesselation
	BRepMesh_IncrementalMesh mesh(shape, deflection);

    std::cout << " -> computing bounds" << std::endl;

	// compute bbox with tesselation
	box.SetVoid();
	BRepBndLib::Add(shape, box);
	box.Get(xMin, yMin, zMin, xMax, yMax, zMax);

	std::cout << " -> done." << std::endl;

    if(argc == 3) {
		std::cout << " -> bounds: " << xMin << ", " << yMin << ", " << zMin << ", " << xMax << ", " << yMax << ", " << zMax << std::endl;
	} else {

		double stlTOL;

		if(argc == 5) {
			stlTOL = parseDouble(argv[4]);
		} else {
			stlTOL = 0.5;
		}

		TopoDS_Shape boundingBox = createBox(xMin,yMin,zMin,xMax,yMax,zMax);
		save(argv[3], boundingBox, stlTOL);

	}
}

// ./occ-csg --transform translate x,y,z         file1.stp file1-translated.stp
void transform(int argc, char *argv[]) {
	if(argc != 6 && argc != 7) {
		error();
	}

	TopoDS_Shape shape = load(argv[4]);

	if(strcmp(argv[2],"translate")==0) {
		std::vector<std::string> values = split(argv[3], ',');

		if(values.size()!=3) {
			error();
		}

		double x1 = parseDouble(values[0]);
		double y1 = parseDouble(values[1]);
		double z1 = parseDouble(values[2]);

		gp_Trsf transformation;
		transformation.SetTranslation(gp_Vec(x1, y1, z1));
		shape = BRepBuilderAPI_GTransform(shape, transformation);
	} else if(strcmp(argv[2],"scale")==0) {
		std::vector<std::string> values = split(argv[3], ',');

		if(values.size()!=3) {
			error();
		}

		double sx = parseDouble(values[0]);
		double sy = parseDouble(values[1]);
		double sz = parseDouble(values[2]);

		gp_GTrsf transformation;

		transformation.SetValue(1,1,sx);
		transformation.SetValue(2,2,sy);
		transformation.SetValue(3,3,sz);
		

		shape = BRepBuilderAPI_GTransform(shape, transformation);
	} else if (strcmp(argv[2],"matrix")==0) {
		std::vector<std::string> values = split(argv[3], ',');

		if(values.size()!=12) {
			error();
		}

		// transformation matrix
		double m[3][4];

		// convert transformation matrix from args
		for(size_t i = 0; i < 3;++i) {
			for(size_t j = 0; j < 4;++j) {
				m[i][j] = parseDouble(values[i+3*j].c_str());
			}
		}

		gp_GTrsf tMat;

		// set transformation matrix
		for(int i = 0; i < 3; i++) {
			for(int j = 0; i < 4; j++) {
				tMat.SetValue(i+1,j+1, m[i][j]);
			}
		}

		shape = BRepBuilderAPI_GTransform(shape, tMat);
	} else {
		error();
	}

	double stlTOL;

	if(argc == 7) {
		stlTOL = parseDouble(argv[6]);
	} else {
		stlTOL = 0.5;
	}

	save(argv[5], shape, stlTOL);
}

TopoDS_Shape createBox(double x1, double y1, double z1, double x2, double y2, double z2) {
	gp_Pnt lowerLeftCornerOfBox(x1,y1,z1);
	gp_Pnt upperRightCornerOfBox(x2,y2,z2);
 	BRepPrimAPI_MakeBox boxMaker(lowerLeftCornerOfBox,upperRightCornerOfBox);
	TopoDS_Shape box = boxMaker.Shape();
	return box;
}

TopoDS_Shape createSphere(double x1, double y1, double z1, double r) {
	gp_Pnt center(x1,y1,z1);
	BRepPrimAPI_MakeSphere sphereMaker(center,r);
	TopoDS_Shape sphere = sphereMaker.Shape();
	return sphere;
}

TopoDS_Shape createCylinder(double r, double h) {
	BRepPrimAPI_MakeCylinder cylinderMaker(r,h);
	TopoDS_Shape cylinder = cylinderMaker.Shape();
	return cylinder;
}

TopoDS_Shape createCylinder(double r, double h, double angle) {
	BRepPrimAPI_MakeCylinder cylinderMaker(r,h,angle);
	TopoDS_Shape cylinder = cylinderMaker.Shape();
	return cylinder;
}

TopoDS_Shape createCone(double r1, double r2, double h) {
	BRepPrimAPI_MakeCone coneMaker(r1,r2,h);
	TopoDS_Shape cone = coneMaker.Shape();
	return cone;
}

TopoDS_Shape createCone(double r1, double r2, double h, double angle) {
	BRepPrimAPI_MakeCone coneMaker(r1,r2,h,angle);
	TopoDS_Shape cone = coneMaker.Shape();
	return cone;
}

TopoDS_Shape extrudePolygon(double ex, double ey, double ez, std::vector<double> const &points) {

    if(points.size()%3!=0) {
		std::cerr << "ERROR: wrong number count, must be multiples of 3, but is " << points.size() << std::endl;
		exit(1);
	}

    size_t numVerts = points.size()/3;
	std::vector<TopoDS_Vertex> vertices;

	for(size_t i = 0; i < points.size(); i+=3) {
		gp_XYZ p;
		p.SetCoord(points[i+0], points[i+1], points[i+2]);

		vertices[i/3] = BRepBuilderAPI_MakeVertex(p);
	}

	BRepBuilderAPI_MakePolygon polyMaker;

	for(size_t i = 0; i < numVerts;i++) {
		polyMaker.Add(vertices[i]);
	}

	polyMaker.Close();

	if(!polyMaker.IsDone()) {
		std::cerr << "ERROR: cannot construct polygon for extrusion. Path invalid (e.g., crossing edges)" << std::endl;
		exit(1);
	}

	TopoDS_Wire wire = polyMaker.Wire();

	if(wire.IsNull()) {
		std::cerr << "ERROR: cannot construct polygon for extrusion. Path invalid (e.g., crossing edges)" << std::endl;
		exit(1);
	}

	TopoDS_Face face = BRepBuilderAPI_MakeFace( wire );

    gp_Vec direction;

	direction.SetX(ex);
	direction.SetY(ey);
	direction.SetZ(ez);

	return BRepPrimAPI_MakePrism(face, direction);
}

TopoDS_Shape extrudeFile(double ex, double ey, double ez, std::string const &filename) {

	TopoDS_Shape face = load(filename);

    gp_Vec direction;

	direction.SetX(ex);
	direction.SetY(ey);
	direction.SetZ(ez);

	return BRepPrimAPI_MakePrism(face, direction);
}

TopoDS_Shape createCircle(double x, double y, double z, double dx, double dy, double dz, double r) {
	gp_Dir dir(dx,dy,dz);
	gp_Pnt point(x,y,z);
	gp_Circ circle(gp_Ax2( point, dir), r);
	BRepBuilderAPI_MakeEdge makeEdge(circle);

    TopoDS_Wire wire = BRepBuilderAPI_MakeWire(makeEdge.Edge());

	TopoDS_Shape shape;

	if( !wire.IsNull()) {
		TopoDS_Shape face = BRepBuilderAPI_MakeFace( wire );
		if(!face.IsNull()) {
			shape = face;
		}
	}

	return shape;
}

TopoDS_Shape createRect2d(double minX, double minY, double maxX, double maxY) {
	std::vector<double> coords;

	coords.push_back(minX); coords.push_back(minY);
	coords.push_back(maxX); coords.push_back(minY);
	coords.push_back(maxX); coords.push_back(maxY);
	coords.push_back(minX); coords.push_back(maxY);

	return createPolygon2d(coords);
}

TopoDS_Shape createPolygon2d(std::vector<double> const &coords) {
	if(coords.size()%2!=0) {
		std::cerr << "ERROR: wrong number count, must be multiples of 2, but is " << coords.size() << std::endl;
		exit(1);
	}

    size_t numVerts = coords.size()/2;
	std::vector<TopoDS_Vertex> vertices;

	for(size_t i = 0; i < coords.size(); i+=2) {
		gp_XYZ p;
		p.SetCoord(coords[i+0], coords[i+1], 0);
		vertices.push_back(BRepBuilderAPI_MakeVertex(p));
	}

	BRepBuilderAPI_MakePolygon polyMaker;

	for(size_t i = 0; i < numVerts;i++) {
		polyMaker.Add(vertices[i]);
	}

	polyMaker.Close();

	if(!polyMaker.IsDone()) {
		std::cerr << "ERROR: cannot construct polygon. Path invalid (e.g., crossing edges)" << std::endl;
		exit(1);
	}

	TopoDS_Wire wire = polyMaker.Wire();

	if(wire.IsNull()) {
		std::cerr << "ERROR: cannot construct polygon. Path invalid (e.g., crossing edges)" << std::endl;
		exit(1);
	}

	TopoDS_Face face = BRepBuilderAPI_MakeFace( wire );

    return face;
}

TopoDS_Shape createText2d(std::string const &font, double fSize, double x, double y, std::string const& text) {

    if(!isAccessible(font)) {
		std::cerr << "ERROR: file '" << font << "' cannot be accessed!" << std::endl;
		exit(1);
	}

	if(!endsWith(toLower(font), ".ttf")) {
		unsupportedFormat(font);
	}

	Font_BRepFont fontObj(font.c_str(), fSize);
	TopoDS_Shape shape = fontObj.RenderText(text.c_str());

	gp_Trsf transformation;
	transformation.SetTranslation(gp_Vec(x, y, 0));
	shape = BRepBuilderAPI_GTransform(shape, transformation);

	return shape;

	//Font_BRepTextBuilder textBuilder;
    //TopoDS_Shape textShape = textBuilder.Perform(font, NCollection_String(text));
}

TopoDS_Shape transform(TopoDS_Shape shape, double transform_matrix[12]) {

	gp_GTrsf tMat;

	// set transformation matrix
	for(int i = 0; i < 3; i++) {
		for(int j = 0; i < 4; j++) {
			tMat.SetValue(i,j, transform_matrix[i+3*j]);
		}
	}

	BRepBuilderAPI_GTransform transform(tMat);

	transform.Perform(shape, true);

	return transform.ModifiedShape(shape);
}

bool isAccessible(std::string const &filename) {
    std::ifstream infile(filename.c_str());
    return infile.good();
}

void unsupportedFormat(std::string const &filename) {
	std::cerr << "> ERROR: unsupported file format: '" << filename << "'" << std::endl;
	exit(1);
}

bool endsWith(std::string const & str, std::string const &ending) {

   if (ending.size() > str.size()) return false;
   return std::equal(ending.rbegin(), ending.rend(), str.rbegin());

}

std::string toLower(std::string const &str) {
	std::string data = str;
	std::transform(data.begin(), data.end(), data.begin(), ::tolower);
	return data;
}

std::vector<std::string> split(std::string const &str, const char sep) {
	std::stringstream ss(str);
	std::vector<std::string> result;

	while( ss.good() )
	{
		std::string substr;
		getline( ss, substr, sep );
		result.push_back( substr );
	}

	return result;
}

double parseDouble(std::string const &str) {
	NUMBER_CONVERSION_ERROR ERR = VALID;
	double res = parseDouble(str, &ERR);

	if(ERR!=VALID) {
		std::cerr << "ERROR: cannot convert number, error_code: " << ERR << std::endl;
		exit(1);
	}

	return res;
}

double parseDouble(std::string const &str, NUMBER_CONVERSION_ERROR *ERROR) {
	const char* buff = str.c_str();
	char *end;
   
   	errno = 0;
   
    const double number = strtod(buff, &end);
   
    if(ERROR != NULL) {
        if (end == buff) {
          *ERROR = ERR_CANNOT_PARSE_NUMBER;
		  return -1;
        } else if ('\0' != *end) {
			*ERROR = ERR_EXTRA_CHARS_AT_THE_END;
          	return -1;
        } else if ((std::numeric_limits<double>::min() == number 
		    || std::numeric_limits<double>::max() == number) && ERANGE == errno) {
            *ERROR = ERR_OUT_OF_RANGE;
			return -1;
        }
    } else {
		std::cerr << "WARNING: performing string-to-double conversion without checks!" << std::endl;
	}

	return number;
}  

int parseInt(std::string const &str) {
	NUMBER_CONVERSION_ERROR ERR = VALID;
	int res = parseInt(str, &ERR);

	if(ERR!=VALID) {
		std::cerr << "ERROR: cannot convert number, error_code: " << ERR << std::endl;
		exit(1);
	}

	return res;
}

int parseInt(std::string const &str, NUMBER_CONVERSION_ERROR *ERROR) {
	const char* buff = str.c_str();
	char *end;
   
   	errno = 0;
   
    const int number = strtol(buff, &end, 10);
   
    if(ERROR != NULL) {
        if (end == buff) {
          *ERROR = ERR_CANNOT_PARSE_NUMBER;
		  return -1;
        } else if ('\0' != *end) {
			*ERROR = ERR_EXTRA_CHARS_AT_THE_END;
          	return -1;
        } else if ((std::numeric_limits<int>::min() == number 
		    || std::numeric_limits<int>::max() == number) && ERANGE == errno) {
            *ERROR = ERR_OUT_OF_RANGE;
			return -1;
        }
    } else {
		std::cerr << "WARNING: performing string-to-int conversion without checks!" << std::endl;
	}

	return number;
}  

void notImplemented() {
	std::cerr << "> ERROR: requested functionality not implemented" << std::endl;
	exit(1);
}

void error() {
	std::cerr << "> ERROR: wrong number of args" << std::endl;
	usage();
	exit(1);
}

void version() {
	std::cout << "Version " << VERSION << std::endl;
}

void usage() {
	std::cerr << "> USAGE: " << std::endl;
	std::cerr << std::endl;
	std::cerr << "Help & Info:" << std::endl;
	std::cerr << std::endl;
	std::cerr << " occ-csg --help" << std::endl;
	std::cerr << " occ-csg --version" << std::endl;
	std::cerr << std::endl;
	std::cerr << "Creating Primitives:" << std::endl;
	std::cerr << std::endl;
	std::cerr << " occ-csg --create box x1,y1,z1,x2,y2,z2                            box.stp" << std::endl;
	std::cerr << " occ-csg --create sphere x1,y1,z1,r                                sphere.stp" << std::endl;
	std::cerr << " occ-csg --create cyl x1,y1,z1,r1,h                                cyl.stp" << std::endl;
	std::cerr << " occ-csg --create 2d:circle x,y,r                                  2dcircle.stp" << std::endl;
	std::cerr << " occ-csg --create 2d:polygon x1,y1,x2,y2,...                       2dpolygon.stp" << std::endl;
	std::cerr << " occ-csg --create 2d:rect x1,y1,x2,y2                              2drectangle.stp" << std::endl;
	std::cerr << " occ-csg --create 2d:text font.ttf 12.0 x,y \"text to render\"       2dtext.stp" << std::endl;
	std::cerr << " occ-csg --create extrusion:polygon ex,ey,ez,x1,y1,z1,x2,y2,z2,... extrude.stp" << std::endl;
	std::cerr << " occ-csg --create extrusion:file ex,ey,ez                          2dpath.stp extrude.stp" << std::endl;
	std::cerr << "" << std::endl;
	std::cerr << "Format Conversion:" << std::endl;
	std::cerr  << std::endl;
	std::cerr << " occ-csg --convert file1.stl file1.stp" << std::endl;
	std::cerr << " occ-csg --convert file1.stp file1.stl 0.1" << std::endl;
	std::cerr << std::endl;
	std::cerr << "Geometric Transformation:" << std::endl;
	std::cerr << std::endl;
	std::cerr << " occ-csg --transform matrix    t1,t2,t3,...,t12 file1.stp file1-transformed.stp" << std::endl;
	std::cerr << " occ-csg --transform translate x,y,z            file1.stp file1-translated.stp" << std::endl;
	std::cerr << " occ-csg --transform scale     sx,sy,sz         file1.stp file1-scaled.stp" << std::endl;
	std::cerr << std::endl;
	std::cerr << "Boolean Operators, Constructive Solid Geometry (CSG):" << std::endl;
	std::cerr << std::endl;
	std::cerr << " occ-csg --csg union file1.stp file2.stp file-out.stp" << std::endl;
	std::cerr << " occ-csg --csg difference file1.stp file2.stp file-out.stp" << std::endl;
	std::cerr << " occ-csg --csg intersection file1.stp file2.stp file-out.stp" << std::endl;
	std::cerr << std::endl;
	std::cerr << "Bounds:" << std::endl;
	std::cerr << std::endl;
	std::cerr << " occ-csg --bounds file.stp" << std::endl;
	std::cerr << " occ-csg --bounds file.stp bounds.stp" << std::endl;
	std::cerr << std::endl;
}
