#include "BackboneRibbonOverlay.h"

#include "Protein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "Atom.h"
#include "Residue.h"
#include "DsspResult.h"
#include "Types.h"

#include <vtkProteinRibbonFilter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkStringArray.h>
#include <vtkIdTypeArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkProperty.h>
#include <vtkNew.h>
#include <vtkTrivialProducer.h>

using namespace nmr;

// Amino acid colors — muted, grouped by physicochemical property.
// Hydrophobic: warm earth tones. Polar: soft blue-green. Charged: muted
// warm (positive) or cool (negative). Aromatic: dusty purple/mauve.
// Special: distinct but not loud.
static void ResidueColor(AminoAcid aa, unsigned char rgb[3]) {
    switch (aa) {
        // Hydrophobic — warm tans, khaki, olive
        case AminoAcid::ALA: rgb[0]=195; rgb[1]=180; rgb[2]=155; break;  // warm sand
        case AminoAcid::VAL: rgb[0]=175; rgb[1]=165; rgb[2]=130; break;  // khaki
        case AminoAcid::LEU: rgb[0]=185; rgb[1]=170; rgb[2]=135; break;  // light khaki
        case AminoAcid::ILE: rgb[0]=165; rgb[1]=155; rgb[2]=120; break;  // olive tan
        case AminoAcid::MET: rgb[0]=190; rgb[1]=180; rgb[2]=110; break;  // muted gold

        // Aromatic — dusty purple/mauve (distinct from everything else)
        case AminoAcid::PHE: rgb[0]=155; rgb[1]=120; rgb[2]=155; break;  // dusty mauve
        case AminoAcid::TRP: rgb[0]=140; rgb[1]=105; rgb[2]=150; break;  // deeper mauve
        case AminoAcid::TYR: rgb[0]=160; rgb[1]=130; rgb[2]=160; break;  // light mauve
        case AminoAcid::HIS: rgb[0]=150; rgb[1]=140; rgb[2]=180; break;  // slate lavender

        // Polar — soft blue-green
        case AminoAcid::SER: rgb[0]=140; rgb[1]=185; rgb[2]=175; break;  // sage
        case AminoAcid::THR: rgb[0]=150; rgb[1]=190; rgb[2]=170; break;  // light sage
        case AminoAcid::ASN: rgb[0]=130; rgb[1]=175; rgb[2]=185; break;  // soft teal
        case AminoAcid::GLN: rgb[0]=140; rgb[1]=180; rgb[2]=190; break;  // light teal
        case AminoAcid::CYS: rgb[0]=190; rgb[1]=190; rgb[2]=120; break;  // muted yellow-green

        // Positive charge — muted warm blue
        case AminoAcid::LYS: rgb[0]=110; rgb[1]=140; rgb[2]=195; break;  // steel blue
        case AminoAcid::ARG: rgb[0]=120; rgb[1]=145; rgb[2]=200; break;  // lighter steel

        // Negative charge — muted terracotta
        case AminoAcid::ASP: rgb[0]=195; rgb[1]=125; rgb[2]=110; break;  // terracotta
        case AminoAcid::GLU: rgb[0]=200; rgb[1]=135; rgb[2]=115; break;  // light terracotta

        // Special
        case AminoAcid::GLY: rgb[0]=210; rgb[1]=210; rgb[2]=210; break;  // light gray
        case AminoAcid::PRO: rgb[0]=190; rgb[1]=160; rgb[2]=145; break;  // dusty rose

        default:             rgb[0]=180; rgb[1]=175; rgb[2]=165; break;  // warm gray
    }
}

BackboneRibbonOverlay::BackboneRibbonOverlay(
    vtkSmartPointer<vtkRenderer> renderer,
    const Protein& protein,
    const ProteinConformation& conf)
    : renderer_(renderer)
{
    // Build the polydata arrays that vtkProteinRibbonFilter expects.
    // Required arrays (from VTK source):
    //   "atom_types"                  vtkStringArray       PDB atom name ("CA","O","N","C",...)
    //   "atom_type"                   vtkIdTypeArray       atomic number
    //   "residue"                     vtkIdTypeArray       residue sequence number
    //   "chain"                       vtkUnsignedCharArray chain ID as char
    //   "secondary_structures"        vtkUnsignedCharArray 'h','s','c' per atom
    //   "secondary_structures_begin"  vtkUnsignedCharArray 1 at SS start
    //   "secondary_structures_end"    vtkUnsignedCharArray 1 at SS end
    //   "ishetatm"                    vtkUnsignedCharArray 0 for protein atoms

    vtkIdType nAtoms = static_cast<vtkIdType>(protein.AtomCount());

    vtkNew<vtkPoints> points;
    points->SetNumberOfPoints(nAtoms);

    vtkNew<vtkStringArray> atomTypes;
    atomTypes->SetName("atom_types");
    atomTypes->SetNumberOfValues(nAtoms);

    vtkNew<vtkIdTypeArray> atomType;
    atomType->SetName("atom_type");
    atomType->SetNumberOfValues(nAtoms);

    vtkNew<vtkIdTypeArray> residue;
    residue->SetName("residue");
    residue->SetNumberOfValues(nAtoms);

    vtkNew<vtkUnsignedCharArray> chain;
    chain->SetName("chain");
    chain->SetNumberOfValues(nAtoms);

    vtkNew<vtkUnsignedCharArray> ss;
    ss->SetName("secondary_structures");
    ss->SetNumberOfValues(nAtoms);

    vtkNew<vtkUnsignedCharArray> ssBegin;
    ssBegin->SetName("secondary_structures_begin");
    ssBegin->SetNumberOfValues(nAtoms);

    vtkNew<vtkUnsignedCharArray> ssEnd;
    ssEnd->SetName("secondary_structures_end");
    ssEnd->SetNumberOfValues(nAtoms);

    vtkNew<vtkUnsignedCharArray> ishetatm;
    ishetatm->SetName("ishetatm");
    ishetatm->SetNumberOfValues(nAtoms);

    bool hasDssp = conf.HasResult<DsspResult>();
    const DsspResult* dssp = hasDssp ? &conf.Result<DsspResult>() : nullptr;

    for (vtkIdType i = 0; i < nAtoms; ++i) {
        size_t si = static_cast<size_t>(i);
        const auto& id = protein.AtomAt(si);
        const auto& res = protein.ResidueAt(id.residue_index);
        Vec3 pos = conf.AtomAt(si).Position();

        points->SetPoint(i, pos.x(), pos.y(), pos.z());
        atomTypes->SetValue(i, id.pdb_atom_name);
        atomType->SetValue(i, static_cast<vtkIdType>(AtomicNumberForElement(id.element)));
        residue->SetValue(i, static_cast<vtkIdType>(res.sequence_number));
        chain->SetValue(i, res.chain_id.empty() ? 'A' : res.chain_id[0]);
        ishetatm->SetValue(i, 0);

        // Secondary structure from DSSP: H=helix, E=sheet, else coil
        unsigned char ssChar = 'c';
        if (dssp) {
            char dsspCode = dssp->SecondaryStructure(id.residue_index);
            if (dsspCode == 'H' || dsspCode == 'G' || dsspCode == 'I')
                ssChar = 'h';
            else if (dsspCode == 'E' || dsspCode == 'B')
                ssChar = 's';
        }
        ss->SetValue(i, ssChar);
        ssBegin->SetValue(i, 0);
        ssEnd->SetValue(i, 0);
    }

    // Mark SS begin/end transitions
    if (dssp) {
        unsigned char prevSS = 'c';
        for (vtkIdType i = 0; i < nAtoms; ++i) {
            unsigned char curSS = ss->GetValue(i);
            if (curSS != prevSS) {
                if (curSS != 'c') ssBegin->SetValue(i, 1);
                if (prevSS != 'c' && i > 0) ssEnd->SetValue(i - 1, 1);
            }
            prevSS = curSS;
        }
        // Last atom ends its SS
        if (nAtoms > 0 && ss->GetValue(nAtoms - 1) != 'c')
            ssEnd->SetValue(nAtoms - 1, 1);
    }

    vtkNew<vtkPolyData> pd;
    pd->SetPoints(points);
    pd->GetPointData()->AddArray(atomTypes);
    pd->GetPointData()->AddArray(atomType);
    pd->GetPointData()->AddArray(residue);
    pd->GetPointData()->AddArray(chain);
    pd->GetPointData()->AddArray(ss);
    pd->GetPointData()->AddArray(ssBegin);
    pd->GetPointData()->AddArray(ssEnd);
    pd->GetPointData()->AddArray(ishetatm);

    // Run the ribbon filter
    vtkNew<vtkTrivialProducer> producer;
    producer->SetOutput(pd);

    vtkNew<vtkProteinRibbonFilter> ribbon;
    ribbon->SetInputConnection(producer->GetOutputPort());
    ribbon->SetCoilWidth(0.5f);
    ribbon->SetHelixWidth(1.8f);
    ribbon->SetSubdivideFactor(20);
    ribbon->SetDrawSmallMoleculesAsSpheres(false);
    ribbon->Update();

    vtkPolyData* ribbonOutput = ribbon->GetOutput();
    if (!ribbonOutput || ribbonOutput->GetNumberOfPoints() == 0)
        return;

    // Replace the filter's default colors with per-residue colors.
    // The filter walks CAs in input order, breaks at chain/residue gaps,
    // and builds strips. Within each strip, each CA contributes subdiv
    // subdivision points, each producing 2 border points (left/right of
    // ribbon) with the same color. We replicate the CA walk to build
    // per-segment residue color lists, then step through the output RGB
    // array in the same order.
    vtkUnsignedCharArray* ribbonRGB = vtkArrayDownCast<vtkUnsignedCharArray>(
        ribbonOutput->GetPointData()->GetScalars("RGB"));
    if (ribbonRGB) {
        int subdiv = ribbon->GetSubdivideFactor();

        // Build per-segment CA lists using the same chain-break logic as
        // vtkProteinRibbonFilter::RequestData (lines 168-189 of VTK source)
        struct Segment { std::vector<AminoAcid> residues; };
        std::vector<Segment> segments;
        segments.push_back({});

        unsigned char prevChainVal = 0;
        vtkIdType prevResiVal = -999;
        bool firstCA = true;

        for (vtkIdType i = 0; i < nAtoms; ++i) {
            std::string atype = atomTypes->GetValue(i);
            if (ishetatm->GetValue(i)) continue;
            if (atype == "CA") {
                unsigned char ch = chain->GetValue(i);
                vtkIdType ri = residue->GetValue(i);
                if (!firstCA && (ch != prevChainVal || ri != prevResiVal + 1)) {
                    // Chain break — start new segment
                    if (segments.back().residues.size() >= 2)
                        segments.push_back({});
                    else
                        segments.back().residues.clear();
                }
                size_t si = static_cast<size_t>(i);
                segments.back().residues.push_back(
                    protein.ResidueAt(protein.AtomAt(si).residue_index).type);
                prevChainVal = ch;
                prevResiVal = ri;
                firstCA = false;
            }
        }

        // Walk the output RGB array, assigning residue colors per segment.
        // Each segment with nCA >= 2 produces a strip with
        // 2 * ((nCA-1)*subdiv + 1) RGB tuples.
        vtkIdType rgbIdx = 0;
        vtkIdType totalTuples = ribbonRGB->GetNumberOfTuples();

        for (const auto& seg : segments) {
            int nCA = static_cast<int>(seg.residues.size());
            if (nCA < 2) continue;
            int len = (nCA - 1) * subdiv + 1;  // subdivision points per strip

            for (int si = 0; si < len; ++si) {
                int caIdx = static_cast<int>(floor(0.5 + si / static_cast<double>(subdiv)));
                if (caIdx >= nCA) caIdx = nCA - 1;

                unsigned char rgb[3];
                ResidueColor(seg.residues[caIdx], rgb);

                // Two border points per subdivision point
                for (int k = 0; k < 2; ++k) {
                    if (rgbIdx < totalTuples)
                        ribbonRGB->SetTypedTuple(rgbIdx++, rgb);
                }
            }
        }
    }

    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputData(ribbonOutput);
    mapper->SetScalarModeToUsePointData();
    mapper->SelectColorArray("RGB");
    mapper->SetColorModeToDirectScalars();

    actor_ = vtkSmartPointer<vtkActor>::New();
    actor_->SetMapper(mapper);
    actor_->GetProperty()->SetInterpolationToPhong();
    actor_->GetProperty()->SetOpacity(0.95);

    renderer_->AddActor(actor_);
}

BackboneRibbonOverlay::~BackboneRibbonOverlay() {
    if (actor_)
        renderer_->RemoveActor(actor_);
}

void BackboneRibbonOverlay::setVisible(bool visible) {
    if (actor_)
        actor_->SetVisibility(visible ? 1 : 0);
}

bool BackboneRibbonOverlay::isVisible() const {
    return actor_ && actor_->GetVisibility();
}
