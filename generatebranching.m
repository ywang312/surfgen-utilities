(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["generatebranching`"];
exportColumbusGeometry::usage="exportColumbusGeometry2[] generates proper COLUMBUS geom file for a set of Cartesian Coordinates.";
generategeom2::usage="generategeom2[] generate Cartesian Coordinates of new geometry based on initial geometry and a displacement vector.";



(* ::Input::Initialization:: *)
exportColumbusGeometry[printdir_,filename_,geometry_List]:=
	Module[{i,cstream,colfile},
		colfile=printdir<>"/"<>filename;
		cstream=OpenWrite[colfile];
               WriteString[cstream,StringJoin[
				" ","N",
				"  ",ToString[PaddedForm[ElementData["Nitrogen","AtomicNumber"],{3,1}]]],
				"  ",ToString[StringForm["``  ``  ``  ``\n",
				         PaddedForm[Part[Part[geometry,1],1],{10,8},ExponentFunction->(Null&)],
					PaddedForm[Part[Part[geometry,1],2],{10,8},ExponentFunction->(Null&)],
					PaddedForm[Part[Part[geometry,1],3],{10,8},ExponentFunction->(Null&)],
					PaddedForm[ElementData["Nitrogen","AtomicWeight"],{10,8}]]]];
               WriteString[cstream,StringJoin[
				" ","C",
				"  ",ToString[PaddedForm[ElementData["Carbon","AtomicNumber"],{3,1}]]],
				"  ",ToString[StringForm["``  ``  ``  ``\n",
				         PaddedForm[Part[Part[geometry,2],1],{10,8},ExponentFunction->(Null&)],
					PaddedForm[Part[Part[geometry,2],2],{10,8},ExponentFunction->(Null&)],
					PaddedForm[Part[Part[geometry,2],3],{10,8},ExponentFunction->(Null&)],
					PaddedForm[ElementData["Carbon","AtomicWeight"],{10,8}]]]];
	       For[i=3,i<=7,i++,
			WriteString[cstream,StringJoin[
                                 " ","H",
				"  ",ToString[PaddedForm[ElementData["Hydrogen","AtomicNumber"],{3,1}]]],
				"  ",ToString[StringForm["``  ``  `` ``\n",
					PaddedForm[Part[Part[geometry,i],1],{10,8},ExponentFunction->(Null&)],
					PaddedForm[Part[Part[geometry,i],2],{10,8},ExponentFunction->(Null&)],
					PaddedForm[Part[Part[geometry,i],3],{10,8},ExponentFunction->(Null&)],
                             PaddedForm[ElementData["Hydrogen","AtomicWeight"],{10,8}]]]]];
		Close[cstream];
];


(* ::Input::Initialization:: *)
generategeom2[initial_,displacement_,numberofgeometry_]:=
             Module[{initialgeom,dispvector,i,inter},
                    initialgeom=Import[initial,"Table"];
                    dispvector=Import[displacement,"Table"];
                    For[i=1,i<=numberofgeometry,i++,
                        inter=initialgeom+dispvector;
                        exportColumbusGeometry[NotebookDirectory[],"geomp"<>ToString[i],inter];
                        initialgeom=inter;]
];


End[]


EndPackage[]
