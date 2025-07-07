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

#include "fuhrmanOmegaWallFunctionFvPatchScalarField.H"
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar fuhrmanOmegaWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void fuhrmanOmegaWallFunctionFvPatchScalarField::checkType()
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

void fuhrmanOmegaWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os)
const
{
    os.writeKeyword("kn") << kn_ << token::END_STATEMENT << nl;
    os.writeKeyword("kr1") << kr1_ << token::END_STATEMENT << nl;
    os.writeKeyword("kr2") << kr2_ << token::END_STATEMENT << nl;
}


void fuhrmanOmegaWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<fuhrmanOmegaWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            fuhrmanOmegaWallFunctionFvPatchScalarField& opf = omegaPatch
            (
                patchi
            );

            if (master == -1)
            {
                master = patchi;
            }

            opf.master() = master;
        }
    }
}


fuhrmanOmegaWallFunctionFvPatchScalarField&
fuhrmanOmegaWallFunctionFvPatchScalarField::omegaPatch(const label patchi)
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const fuhrmanOmegaWallFunctionFvPatchScalarField& opf =
        refCast<const fuhrmanOmegaWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<fuhrmanOmegaWallFunctionFvPatchScalarField&>(opf);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fuhrmanOmegaWallFunctionFvPatchScalarField::
fuhrmanOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    kn_(1e-6),
    kr1_(200e0),
    kr2_(100e0),
    omega_(),
    initialised_(false),
    master_(-1)
{
    checkType();
}


fuhrmanOmegaWallFunctionFvPatchScalarField::
fuhrmanOmegaWallFunctionFvPatchScalarField
(
    const fuhrmanOmegaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    kn_(ptf.kn_),
    kr1_(ptf.kr1_),
    kr2_(ptf.kr2_),
    omega_(),
    initialised_(false),
    master_(-1)
{
    checkType();
}


fuhrmanOmegaWallFunctionFvPatchScalarField::
fuhrmanOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    kn_(dict.lookupOrDefault<scalar>("kn", 1e-6)),
    kr1_(dict.lookupOrDefault<scalar>("kr1", 200e0)),
    kr2_(dict.lookupOrDefault<scalar>("kr2", 100e0)),
    omega_(),
    initialised_(false),
    master_(-1)
{
    checkType();

    // apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


fuhrmanOmegaWallFunctionFvPatchScalarField::
fuhrmanOmegaWallFunctionFvPatchScalarField
(
    const fuhrmanOmegaWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedValueFvPatchField<scalar>(owfpsf),
    kn_(owfpsf.kn_),
    kr1_(owfpsf.kr1_),
    kr2_(owfpsf.kr2_),
    omega_(),
    initialised_(false),
    master_(-1)
{
    checkType();
}


fuhrmanOmegaWallFunctionFvPatchScalarField::
fuhrmanOmegaWallFunctionFvPatchScalarField
(
    const fuhrmanOmegaWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    kn_(owfpsf.kn_),
    kr1_(owfpsf.kr1_),
    kr2_(owfpsf.kr2_),
    omega_(),
    initialised_(false),
    master_(-1)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalarField& fuhrmanOmegaWallFunctionFvPatchScalarField::omega(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            omega_ = 0.0;
        }

        return omega_;
    }

    return omegaPatch(master_).omega(init);
}


void fuhrmanOmegaWallFunctionFvPatchScalarField::updateCoeffs()
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

    setMaster();

    scalarField newOmegaValues(this->patchInternalField());

    const label patchi = patch().index();

    const tmp<volScalarField> tk = turbModel.k();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const tmp<volScalarField> tnut = turbModel.nut();
    const scalarField& nutw = tnutw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    // Set omega
    forAll(nutw, facei)
    {
        scalar utau = max(1.0e-6, sqrt((nutw[facei] + nuw[facei])*
                      magGradUw[facei]));
        scalar knplus = kn_*utau/nuw[facei];
        scalar SR = 0e0;
        if (knplus<=5e0)
        {
            SR = pow(kr1_/knplus, 2);
        }
        else
        {
            SR = kr2_/knplus + (pow(kr1_/knplus, 2)-kr2_/knplus)*exp(5-knplus);
        }

        newOmegaValues[facei] = pow(utau, 2)*SR/nuw[facei];
    }
    this->operator==(newOmegaValues);

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void fuhrmanOmegaWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeLocalEntries(os);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fuhrmanOmegaWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
