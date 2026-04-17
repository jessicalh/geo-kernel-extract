#include "QtBackboneRibbonOverlay.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/QtFrame.h"
#include "../model/QtResidue.h"

#include <QLoggingCategory>

#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkTrivialProducer.h>

#include <cmath>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cRibbon, "h5reader.overlay.ribbon")

// Per-amino-acid ribbon colour. Typed dispatch on AminoAcid enum —
// never a string comparison. Palette matches the existing ui/ viewer
// (muted, grouped by physicochemical property).
void ResidueRgb(model::AminoAcid aa, unsigned char rgb[3]) {
    switch (aa) {
        // Hydrophobic — warm tans/khaki/olive
        case model::AminoAcid::ALA: rgb[0]=195; rgb[1]=180; rgb[2]=155; break;
        case model::AminoAcid::VAL: rgb[0]=175; rgb[1]=165; rgb[2]=130; break;
        case model::AminoAcid::LEU: rgb[0]=185; rgb[1]=170; rgb[2]=135; break;
        case model::AminoAcid::ILE: rgb[0]=165; rgb[1]=155; rgb[2]=120; break;
        case model::AminoAcid::MET: rgb[0]=190; rgb[1]=180; rgb[2]=110; break;
        // Aromatic — dusty mauve
        case model::AminoAcid::PHE: rgb[0]=155; rgb[1]=120; rgb[2]=155; break;
        case model::AminoAcid::TRP: rgb[0]=140; rgb[1]=105; rgb[2]=150; break;
        case model::AminoAcid::TYR: rgb[0]=160; rgb[1]=130; rgb[2]=160; break;
        case model::AminoAcid::HIS: rgb[0]=150; rgb[1]=140; rgb[2]=180; break;
        // Polar — soft blue-green
        case model::AminoAcid::SER: rgb[0]=140; rgb[1]=185; rgb[2]=175; break;
        case model::AminoAcid::THR: rgb[0]=150; rgb[1]=190; rgb[2]=170; break;
        case model::AminoAcid::ASN: rgb[0]=130; rgb[1]=175; rgb[2]=185; break;
        case model::AminoAcid::GLN: rgb[0]=140; rgb[1]=180; rgb[2]=190; break;
        case model::AminoAcid::CYS: rgb[0]=190; rgb[1]=190; rgb[2]=120; break;
        // Positive — steel blue
        case model::AminoAcid::LYS: rgb[0]=110; rgb[1]=140; rgb[2]=195; break;
        case model::AminoAcid::ARG: rgb[0]=120; rgb[1]=145; rgb[2]=200; break;
        // Negative — terracotta
        case model::AminoAcid::ASP: rgb[0]=195; rgb[1]=125; rgb[2]=110; break;
        case model::AminoAcid::GLU: rgb[0]=200; rgb[1]=135; rgb[2]=115; break;
        // Special
        case model::AminoAcid::GLY: rgb[0]=210; rgb[1]=210; rgb[2]=210; break;
        case model::AminoAcid::PRO: rgb[0]=190; rgb[1]=160; rgb[2]=145; break;
        default:                     rgb[0]=180; rgb[1]=175; rgb[2]=165; break;
    }
}

// Map DsspCode (8-class) to vtkProteinRibbonFilter's input 'h' (helix),
// 's' (sheet), or 'c' (coil).
unsigned char SsCharForDssp(model::DsspCode c) {
    switch (c) {
        case model::DsspCode::AlphaHelix:
        case model::DsspCode::Helix310:
        case model::DsspCode::PiHelix:   return 'h';
        case model::DsspCode::BetaStrand:
        case model::DsspCode::BetaBridge: return 's';
        default:                          return 'c';
    }
}

}  // namespace

QtBackboneRibbonOverlay::QtBackboneRibbonOverlay(
    vtkSmartPointer<vtkRenderer> renderer, QObject* parent)
    : QObject(parent), renderer_(std::move(renderer))
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("QtBackboneRibbonOverlay"));
}

QtBackboneRibbonOverlay::~QtBackboneRibbonOverlay() {
    if (actor_) renderer_->RemoveActor(actor_);
}

void QtBackboneRibbonOverlay::Build(
    const model::QtProtein&      protein,
    const model::QtConformation& conformation)
{
    ASSERT_THREAD(this);

    if (protein_ == &protein && conformation_ == &conformation && actor_)
        return;

    if (actor_) renderer_->RemoveActor(actor_);
    actor_ = nullptr;

    protein_      = &protein;
    conformation_ = &conformation;

    const vtkIdType nAtoms = static_cast<vtkIdType>(protein.atomCount());

    points_ = vtkSmartPointer<vtkPoints>::New();
    points_->SetNumberOfPoints(nAtoms);

    atomTypes_ = vtkSmartPointer<vtkStringArray>::New();
    atomTypes_->SetName("atom_types");
    atomTypes_->SetNumberOfValues(nAtoms);

    atomType_ = vtkSmartPointer<vtkIdTypeArray>::New();
    atomType_->SetName("atom_type");
    atomType_->SetNumberOfValues(nAtoms);

    residue_ = vtkSmartPointer<vtkIdTypeArray>::New();
    residue_->SetName("residue");
    residue_->SetNumberOfValues(nAtoms);

    chain_ = vtkSmartPointer<vtkUnsignedCharArray>::New();
    chain_->SetName("chain");
    chain_->SetNumberOfValues(nAtoms);

    ss_ = vtkSmartPointer<vtkUnsignedCharArray>::New();
    ss_->SetName("secondary_structures");
    ss_->SetNumberOfValues(nAtoms);

    ssBegin_ = vtkSmartPointer<vtkUnsignedCharArray>::New();
    ssBegin_->SetName("secondary_structures_begin");
    ssBegin_->SetNumberOfValues(nAtoms);

    ssEnd_ = vtkSmartPointer<vtkUnsignedCharArray>::New();
    ssEnd_->SetName("secondary_structures_end");
    ssEnd_->SetNumberOfValues(nAtoms);

    ishetatm_ = vtkSmartPointer<vtkUnsignedCharArray>::New();
    ishetatm_->SetName("ishetatm");
    ishetatm_->SetNumberOfValues(nAtoms);

    // Identity-only columns (constant across the trajectory).
    for (vtkIdType i = 0; i < nAtoms; ++i) {
        const size_t si = static_cast<size_t>(i);
        const auto& atom = protein.atom(si);
        const auto& res  = protein.residue(atom.residueIndex);

        atomTypes_->SetValue(i, atom.h5AtomName.toStdString());
        atomType_->SetValue(i,
            static_cast<vtkIdType>(model::AtomicNumberForElement(atom.element)));
        residue_->SetValue(i, static_cast<vtkIdType>(res.residueNumber));
        chain_->SetValue(i,
            res.chainId.isEmpty() ? 'A' : res.chainId.at(0).toLatin1());
        ishetatm_->SetValue(i, 0);
    }

    inputPd_ = vtkSmartPointer<vtkPolyData>::New();
    inputPd_->SetPoints(points_);
    inputPd_->GetPointData()->AddArray(atomTypes_);
    inputPd_->GetPointData()->AddArray(atomType_);
    inputPd_->GetPointData()->AddArray(residue_);
    inputPd_->GetPointData()->AddArray(chain_);
    inputPd_->GetPointData()->AddArray(ss_);
    inputPd_->GetPointData()->AddArray(ssBegin_);
    inputPd_->GetPointData()->AddArray(ssEnd_);
    inputPd_->GetPointData()->AddArray(ishetatm_);

    auto producer = vtkSmartPointer<vtkTrivialProducer>::New();
    producer->SetOutput(inputPd_);

    ribbon_ = vtkSmartPointer<vtkProteinRibbonFilter>::New();
    ribbon_->SetInputConnection(producer->GetOutputPort());
    ribbon_->SetCoilWidth(0.5f);
    ribbon_->SetHelixWidth(1.8f);
    ribbon_->SetSubdivideFactor(subdivideFactor_);
    ribbon_->SetDrawSmallMoleculesAsSpheres(false);

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(ribbon_->GetOutputPort());
    mapper->SetScalarModeToUsePointData();
    mapper->SelectColorArray("RGB");
    mapper->SetColorModeToDirectScalars();

    actor_ = vtkSmartPointer<vtkActor>::New();
    actor_->SetMapper(mapper);
    actor_->GetProperty()->SetInterpolationToPhong();
    actor_->GetProperty()->SetOpacity(0.95);

    renderer_->AddActor(actor_);

    qCInfo(cRibbon).noquote()
        << "Built ribbon pipeline | atoms=" << nAtoms
        << "| subdivide=" << subdivideFactor_;

    setFrame(0);
    setVisible(visible_);
}

void QtBackboneRibbonOverlay::UpdateInputArrays(int t) {
    const size_t tIndex = static_cast<size_t>(t);
    const auto& frame = conformation_->frame(tIndex);
    const vtkIdType nAtoms = static_cast<vtkIdType>(protein_->atomCount());

    for (vtkIdType i = 0; i < nAtoms; ++i) {
        const size_t si = static_cast<size_t>(i);
        const model::Vec3 p = frame.position(si);
        points_->SetPoint(i, p.x(), p.y(), p.z());

        const auto& atom = protein_->atom(si);
        const model::DsspCode code =
            atom.residueIndex >= 0
                ? frame.dsspCode(static_cast<size_t>(atom.residueIndex))
                : model::DsspCode::Unknown;
        ss_->SetValue(i, SsCharForDssp(code));
        ssBegin_->SetValue(i, 0);
        ssEnd_->SetValue(i, 0);
    }

    // Mark SS begin/end transitions by scanning the per-atom SS string.
    unsigned char prevSS = 'c';
    for (vtkIdType i = 0; i < nAtoms; ++i) {
        const unsigned char curSS = ss_->GetValue(i);
        if (curSS != prevSS) {
            if (curSS != 'c')                ssBegin_->SetValue(i, 1);
            if (prevSS != 'c' && i > 0)      ssEnd_->SetValue(i - 1, 1);
        }
        prevSS = curSS;
    }
    if (nAtoms > 0 && ss_->GetValue(nAtoms - 1) != 'c')
        ssEnd_->SetValue(nAtoms - 1, 1);

    points_->Modified();
    ss_->Modified();
    ssBegin_->Modified();
    ssEnd_->Modified();
    inputPd_->Modified();
}

void QtBackboneRibbonOverlay::ApplyResidueColors() {
    vtkPolyData* ribbonOutput = ribbon_->GetOutput();
    if (!ribbonOutput || ribbonOutput->GetNumberOfPoints() == 0) return;

    auto* rgbArray = vtkArrayDownCast<vtkUnsignedCharArray>(
        ribbonOutput->GetPointData()->GetScalars("RGB"));
    if (!rgbArray) return;

    // Rebuild per-segment residue lists using the same chain-break logic
    // as vtkProteinRibbonFilter::RequestData (CA atoms only, break on
    // chain change or non-consecutive residue number).
    segments_.clear();
    segments_.push_back({});

    unsigned char prevChainVal = 0;
    vtkIdType     prevResiVal  = -999;
    bool          firstCA      = true;

    const vtkIdType nAtoms = static_cast<vtkIdType>(protein_->atomCount());
    for (vtkIdType i = 0; i < nAtoms; ++i) {
        if (ishetatm_->GetValue(i)) continue;
        const std::string atype = atomTypes_->GetValue(i);
        // CA dispatch is by NAME here because vtkProteinRibbonFilter
        // itself walks CAs by name — we must mirror its traversal to
        // align our segment boundaries with its output point layout.
        // This is the ribbon filter's contract; we are staying inside
        // its lines, not dispatching physics on names elsewhere.
        if (atype != "CA") continue;

        const unsigned char ch = chain_->GetValue(i);
        const vtkIdType     ri = residue_->GetValue(i);
        if (!firstCA && (ch != prevChainVal || ri != prevResiVal + 1)) {
            if (segments_.back().residues.size() >= 2)
                segments_.push_back({});
            else
                segments_.back().residues.clear();
        }
        const size_t si = static_cast<size_t>(i);
        const auto& res = protein_->residue(protein_->atom(si).residueIndex);
        segments_.back().residues.push_back(res.aminoAcid);
        prevChainVal = ch;
        prevResiVal  = ri;
        firstCA      = false;
    }

    vtkIdType rgbIdx     = 0;
    const vtkIdType total = rgbArray->GetNumberOfTuples();
    for (const auto& seg : segments_) {
        const int nCA = static_cast<int>(seg.residues.size());
        if (nCA < 2) continue;
        const int len = (nCA - 1) * subdivideFactor_ + 1;
        for (int si = 0; si < len; ++si) {
            int caIdx = static_cast<int>(
                std::floor(0.5 + si / static_cast<double>(subdivideFactor_)));
            if (caIdx >= nCA) caIdx = nCA - 1;
            unsigned char rgb[3];
            ResidueRgb(seg.residues[caIdx], rgb);
            for (int k = 0; k < 2; ++k) {
                if (rgbIdx < total) rgbArray->SetTypedTuple(rgbIdx++, rgb);
            }
        }
    }
    rgbArray->Modified();
}

void QtBackboneRibbonOverlay::setFrame(int t) {
    ASSERT_THREAD(this);
    if (!protein_ || !conformation_) return;
    if (t < 0 || static_cast<size_t>(t) >= conformation_->frameCount()) return;

    UpdateInputArrays(t);
    ribbon_->Update();
    ApplyResidueColors();
    // MoleculeScene renders after all overlays update this frame.
}

void QtBackboneRibbonOverlay::setVisible(bool visible) {
    ASSERT_THREAD(this);
    visible_ = visible;
    if (actor_) actor_->SetVisibility(visible ? 1 : 0);
}

}  // namespace h5reader::app
