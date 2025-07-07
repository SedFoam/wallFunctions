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

#include "leeOmegaWallFunctionFvPatchScalarField.H"
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

void leeOmegaWallFunctionFvPatchScalarField::checkType()
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

void leeOmegaWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os)
const
{
    os.writeKeyword("kn") << kn_ << token::END_STATEMENT << nl;
    os.writeEntryIfDifferent<scalar>("beta1", 0.075, beta1_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

leeOmegaWallFunctionFvPatchScalarField::
leeOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    kn_(1e-6),
    d0_(0.03*kn_),
    beta1_(0.075),
    wallCoeffs_()
{
    checkType();
}


leeOmegaWallFunctionFvPatchScalarField::
leeOmegaWallFunctionFvPatchScalarField
(
    const leeOmegaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    kn_(ptf.kn_),
    d0_(0.03*kn_),
    beta1_(ptf.beta1_),
    wallCoeffs_(ptf.wallCoeffs_)
{
    checkType();
}


leeOmegaWallFunctionFvPatchScalarField::
leeOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    kn_(dict.lookupOrDefault<scalar>("kn", 1e-6)),
    d0_(0.03*kn_),
    beta1_(dict.getOrDefault<scalar>("beta1", 0.075)),
    wallCoeffs_(dict)
{
    checkType();

    // apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


leeOmegaWallFunctionFvPatchScalarField::
leeOmegaWallFunctionFvPatchScalarField
(
    const leeOmegaWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedValueFvPatchField<scalar>(owfpsf),
    kn_(owfpsf.kn_),
    d0_(0.03*kn_),
    beta1_(owfpsf.beta1_),
    wallCoeffs_(owfpsf.wallCoeffs_)
{
    checkType();
}


leeOmegaWallFunctionFvPatchScalarField::
leeOmegaWallFunctionFvPatchScalarField
(
    const leeOmegaWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    kn_(owfpsf.kn_),
    d0_(0.03*kn_),
    beta1_(owfpsf.beta1_),
    wallCoeffs_(owfpsf.wallCoeffs_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void leeOmegaWallFunctionFvPatchScalarField::updateCoeffs()
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

    scalarField newOmegaValues(this->patchInternalField());

    const label patchi = patch().index();

    // access betaStar constant from turbulence model
    const scalar sqrtBetaStar = sqrt(wallCoeffs_.Cmu());
    // access von Karman constant from turbulence model
    const scalar kappa = wallCoeffs_.kappa();

    // first cell center distance to wall patch
    const scalarField& y = turbModel.y()[patchi];

    // fluid kinematic viscosity
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const tmp<volScalarField> tnut = turbModel.nut();
    const scalarField& nutw = tnutw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    scalar smallVal = SMALL;

    // Set omega
    forAll(nutw, facei)
    {
        // friction velocity
        scalar utau = max(1.0e-6, sqrt((nutw[facei] + nuw[facei])*
                      magGradUw[facei]));
        // roughness Reynolds number
        scalar knPlus = kn_*utau/nuw[facei];

        // modified hydrodynamical roughness
        scalar d0Tilde =
            d0_ * min(1, pow(knPlus/30, 2./3.))
            * min(1, pow(knPlus/45, 0.25))
            * min(1, pow(knPlus/60, 0.25));

        // omega value on face, Lee (2010)
        newOmegaValues[facei] = min(
            (utau * Foam::log(1 + y[facei]/d0Tilde)) / (
                sqrtBetaStar * kappa * (d0Tilde + smallVal)),
            6 * nuw[facei] / (beta1_ * sqr(y[facei]))
        );
    }
    this->operator==(newOmegaValues);

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void leeOmegaWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeLocalEntries(os);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    leeOmegaWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
