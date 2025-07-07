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

\*---------------------------------------------------------------------------*/

#include "knoppkWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "nutkWallFunctionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void knoppkWallFunctionFvPatchScalarField::checkType()
{
    if (not isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}

void knoppkWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os)
const
{
    os.writeKeyword("kn") << kn_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

knoppkWallFunctionFvPatchScalarField::
knoppkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    kn_(1e-6),
    wallCoeffs_(),
    initialised_(false)
{
    checkType();
}


knoppkWallFunctionFvPatchScalarField::
knoppkWallFunctionFvPatchScalarField
(
    const knoppkWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    kn_(ptf.kn_),
    wallCoeffs_(ptf.wallCoeffs_),
    initialised_(false)
{
    checkType();
}


knoppkWallFunctionFvPatchScalarField::
knoppkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    kn_(dict.lookupOrDefault<scalar>("kn", 1e-6)),
    wallCoeffs_(dict),
    initialised_(false)
{
    checkType();

    // apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


knoppkWallFunctionFvPatchScalarField::
knoppkWallFunctionFvPatchScalarField
(
    const knoppkWallFunctionFvPatchScalarField& kwfpsf
)
:
    fixedValueFvPatchField<scalar>(kwfpsf),
    kn_(kwfpsf.kn_),
    wallCoeffs_(kwfpsf.wallCoeffs_),
    initialised_(false)
{
    checkType();
}


knoppkWallFunctionFvPatchScalarField::
knoppkWallFunctionFvPatchScalarField
(
    const knoppkWallFunctionFvPatchScalarField& kwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(kwfpsf, iF),
    kn_(kwfpsf.kn_),
    wallCoeffs_(kwfpsf.wallCoeffs_),
    initialised_(false)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void knoppkWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    scalarField newkValues(this->patchInternalField());

    const label patchi = patch().index();

    // access betaStar coefficient from turbulence model
    const scalar sqrtBetaStar = sqrt(wallCoeffs_.Cmu());

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const tmp<volScalarField> tnut = turbModel.nut();
    const scalarField& nutw = tnutw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    // Set k values on patch
    forAll(nutw, facei)
    {
        // friction velocity
        scalar utau = max(1.0e-6, sqrt((nutw[facei] + nuw[facei])*
                      magGradUw[facei]));
        // rugh Reynolds number
        scalar knPlus = kn_*utau/nuw[facei];

        // new value for k on patch face, knopp (2010)
        newkValues[facei] =
            min(1., knPlus / 90.)
            * sqr(utau) / sqrtBetaStar;
    }
    this->operator==(newkValues);

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void knoppkWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeLocalEntries(os);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    knoppkWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
