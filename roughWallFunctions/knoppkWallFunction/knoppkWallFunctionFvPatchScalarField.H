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
    Foam::knoppkWallFunctionFvPatchScalarField

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
    knoppkWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef knoppkWallFunctionFvPatchScalarField_H
#define knoppkWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"
#include "wallFunctionCoefficients.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class turbulenceModel;

/*---------------------------------------------------------------------------*\
             Class knoppkWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class knoppkWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
protected:

    // Protected data

    //- Tolerance used in weighted calculations
    static scalar tolerance_;

    //- kn coefficient (roughness)
    scalar kn_;

    //- Wall-function coefficients
    wallFunctionCoefficients wallCoeffs_;

    //- Initialised flag
    bool initialised_;


    // Protected Member Functions

        //- Check the type of the patch
        virtual void checkType();

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("knoppkWallFunction");


    // Constructors

    //- Construct from patch and internal field
    knoppkWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    knoppkWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );
    
    //- Construct by mapping given
    // knoppkWallFunctionFvPatchScalarField
    //  onto a new patch
    knoppkWallFunctionFvPatchScalarField
    (
        const knoppkWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    knoppkWallFunctionFvPatchScalarField
    (
        const knoppkWallFunctionFvPatchScalarField&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
            (
                new knoppkWallFunctionFvPatchScalarField(*this)
            );
    }

    //- Construct as copy setting internal field reference
    knoppkWallFunctionFvPatchScalarField
    (
        const knoppkWallFunctionFvPatchScalarField&,
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
                new knoppkWallFunctionFvPatchScalarField(*this, iF)
            );
    }


    // Member functions

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
