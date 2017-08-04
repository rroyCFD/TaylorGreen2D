/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    TaylorGreen2D

Description
    Initializes the flow field for 2D Taylor Green Vortex.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readTransportProperties.H"
    
    // Read in the existing solution files.   
    Info << "Reading field U" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    
    Info << "Reading field p" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    
    dimensionedScalar Re = Uinit*L/nu;
    Info<< "\nInitializing Taylor-Greem Vortex flow with Re = "<< Re.value() << "\n" << endl;
    
    // Intialize Velocity field
    volScalarField x = mesh.C().component(vector::X);
    volScalarField y = mesh.C().component(vector::Y);
    volScalarField z = mesh.C().component(vector::Z); 
    
    U  =  Uinit* (vector(1,0,0) * sin(x/L) * cos(y/L) * cos(z/L) 
                - vector(0,1,0) * cos(x/L) * sin(y/L) * cos(z/L) 
                + vector(0,0,1) * scalar(0.));
    U.correctBoundaryConditions();
    Info<< "Writing field U" << endl;
    U.write();
    
    p = sqr(Uinit)/16. * (cos(2.*x/L) + cos(2.*y/L))*(cos(2.*z/L)+scalar(2.));
    p.correctBoundaryConditions();
    Info<< "Writing field p" << endl;
    p.write();
 
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
