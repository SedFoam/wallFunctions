/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::fuhrmanOmegaWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall function constraint on turbulnce
    specific dissipation, omega.  The values are computed using:

        \f[
            \omega = sqrt(\omega_{vis}^2 + \omega_{log}^2)
        \f]

    where

    \vartable
        \omega_{vis} | omega in viscous region
        \omega_{log} | omega in logarithmic region
    \endvartable

    Model described by Eq.(15) of:
    \verbatim
        Menter, F., Esch, T.
        "Elements of Industrial Heat Transfer Prediction"
        16th Brazilian Congress of Mechanical Engineering (COBEM),
        Nov. 2001
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        Cmu          | model coefficient       | no          | 0.09
        kappa        | Von Karman constant     | no          | 0.41
        E            | model coefficient       | no          | 9.8
        beta1        | model coefficient       | no          | 0.075
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            omegaWallFunction;
    }
    \endverbatim

SourceFiles
    fuhrmanOmegaWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fuhrmanOmegaWallFunctionFvPatchScalarField_H
#define fuhrmanOmegaWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class turbulenceModel;

/*---------------------------------------------------------------------------*\
             Class fuhrmanOmegaWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fuhrmanOmegaWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
protected:

    // Protected data

        //- Tolerance used in weighted calculations
        static scalar tolerance_;

         //- kn coefficient (roughness)
        scalar kn_;

        //- kr1 empirical coefficient (smooth)
        scalar kr1_;

        //- kr2 empirical coefficient (rough)
        scalar kr2_;

        //- Local copy of turbulence omega field
        scalarField omega_;

        //- Initialised flag
        bool initialised_;

        //- Master patch ID
        label master_;


    // Protected Member Functions

        //- Check the type of the patch
        virtual void checkType();

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;

        //- Set the master patch - master is responsible for updating all
        //  wall function patches
        virtual void setMaster();

        //- Helper function to return non-const access to an omega patch
        virtual fuhrmanOmegaWallFunctionFvPatchScalarField& omegaPatch
        (
            const label patchi
        );

        //- Return non-const access to the master patch ID
        virtual label& master()
        {
            return master_;
        }


public:

    //- Runtime type information
    TypeName("fuhrmanOmegaWallFunction");


    // Constructors

        //- Construct from patch and internal field
        fuhrmanOmegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fuhrmanOmegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // fuhrmanOmegaWallFunctionFvPatchScalarField
        //  onto a new patch
        fuhrmanOmegaWallFunctionFvPatchScalarField
        (
            const fuhrmanOmegaWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fuhrmanOmegaWallFunctionFvPatchScalarField
        (
            const fuhrmanOmegaWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fuhrmanOmegaWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fuhrmanOmegaWallFunctionFvPatchScalarField
        (
            const fuhrmanOmegaWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new fuhrmanOmegaWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return non-const access to the master's omega field
            scalarField& omega(bool init = false);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
