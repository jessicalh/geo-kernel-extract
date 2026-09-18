#include "LearnedActivityDock.h"

#include "MoleculeScene.h"
#include "QtPlaybackController.h"
#include "TensorGlyphActor.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/QtProtein.h"
#include "../model/TransformedConformation.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFile>
#include <QFileDialog>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QJsonDocument>
#include <QLabel>
#include <QLoggingCategory>
#include <QMessageBox>
#include <QSignalBlocker>
#include <QStyle>
#include <QToolButton>
#include <QVBoxLayout>

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace h5reader::app {
namespace {
Q_LOGGING_CATEGORY(cActivity, "h5reader.learned_activity")
constexpr double kPositionToleranceA = 0.001;
constexpr double kTraceTolerance = 1e-5;

void require(bool condition, const QString& message) {
    if (!condition)
        throw std::runtime_error(message.toStdString());
}

int indexValue(const QJsonValue& value, const QString& name) {
    const double number = value.toDouble(-1.0);
    require(value.isDouble() && std::isfinite(number) && number >= 0.0
                && number <= std::numeric_limits<int>::max() && std::floor(number) == number,
            name + QStringLiteral(" must be a nonnegative integer"));
    return static_cast<int>(number);
}

std::vector<double> numbers(const QJsonValue& value, int count) {
    const auto array = value.toArray();
    require(value.isArray() && array.size() == count, QStringLiteral("wrong numeric array width"));
    std::vector<double> result;
    for (const auto& entry : array) {
        require(entry.isDouble() && std::isfinite(entry.toDouble()), QStringLiteral("nonfinite or nonnumeric value"));
        result.push_back(entry.toDouble());
    }
    return result;
}

QJsonArray vectorJson(const model::Vec3& value) {
    return {value.x(), value.y(), value.z()};
}

QJsonArray tensorJson(const model::Mat3& value) {
    return {value(0, 0), value(0, 1), value(0, 2), value(1, 1), value(1, 2), value(2, 2)};
}
}  // namespace

LearnedActivityDock::LearnedActivityDock(QWidget* parent)
    : QDockWidget(QStringLiteral("Learned activity"), parent) {
    setObjectName(QStringLiteral("LearnedActivityDock"));
    auto* body = new QWidget(this);
    auto* layout = new QVBoxLayout(body);
    auto* buttons = new QHBoxLayout;
    auto* open = new QToolButton(body);
    open->setIcon(style()->standardIcon(QStyle::SP_DialogOpenButton));
    open->setToolTip(QStringLiteral("Open activation capture"));
    open->setAccessibleName(open->toolTip());
    auto* remove = new QToolButton(body);
    remove->setIcon(style()->standardIcon(QStyle::SP_DialogResetButton));
    remove->setToolTip(QStringLiteral("Clear activation capture"));
    remove->setAccessibleName(remove->toolTip());
    visible_ = new QCheckBox(QStringLiteral("Show surfaces"), body);
    visible_->setChecked(true);
    buttons->addWidget(open);
    buttons->addWidget(remove);
    buttons->addWidget(visible_);
    buttons->addStretch();
    layout->addLayout(buttons);
    auto* form = new QFormLayout;
    channel_ = new QComboBox(body);
    channel_->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    channel_->setMinimumContentsLength(20);
    radius_ = new QDoubleSpinBox(body);
    radius_->setRange(0.01, 5.0);
    radius_->setDecimals(2);
    radius_->setSingleStep(0.1);
    radius_->setValue(1.5);
    radius_->setToolTip(QStringLiteral("Angstroms at the largest absolute tensor eigenvalue in this entire capture. Shared across atoms, channels and frames."));
    opacity_ = new QDoubleSpinBox(body);
    opacity_->setRange(0.0, 1.0);
    opacity_->setSingleStep(0.05);
    opacity_->setValue(0.65);
    form->addRow(QStringLiteral("Channel"), channel_);
    form->addRow(QStringLiteral("Peak radius (A)"), radius_);
    form->addRow(QStringLiteral("Opacity"), opacity_);
    layout->addLayout(form);
    status_ = new QLabel(body);
    status_->setWordWrap(true);
    status_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(status_);
    auto* key = new QLabel(QStringLiteral("Red: negative   Blue: positive\nHidden 2e activity; model units"), body);
    key->setWordWrap(true);
    layout->addWidget(key);
    layout->addStretch();
    setWidget(body);
    connect(open, &QToolButton::clicked, this, &LearnedActivityDock::openCapture);
    connect(remove, &QToolButton::clicked, this, &LearnedActivityDock::clear);
    connect(channel_, &QComboBox::currentIndexChanged, this, &LearnedActivityDock::refresh);
    connect(radius_, &QDoubleSpinBox::valueChanged, this, &LearnedActivityDock::refresh);
    connect(opacity_, &QDoubleSpinBox::valueChanged, this, &LearnedActivityDock::refresh);
    connect(visible_, &QCheckBox::toggled, this, &LearnedActivityDock::refresh);
    updateControls();
}

LearnedActivityDock::~LearnedActivityDock() = default;

void LearnedActivityDock::setContext(MoleculeScene* scene, const model::QtProtein* protein,
                                     model::Conformation* raw, model::TransformedConformation* display,
                                     QtPlaybackController* playback) {
    ASSERT_THREAD(this);
    disconnect(frameConnection_);
    disconnect(transformConnection_);
    clear();
    scene_ = scene;
    protein_ = protein;
    raw_ = raw;
    display_ = display;
    playback_ = playback;
    if (playback_)
        frameConnection_ = connect(playback, &QtPlaybackController::frameChanged, this, &LearnedActivityDock::refresh);
    if (display_)
        transformConnection_ = connect(display, &model::TransformedConformation::transformChanged,
                                       this, &LearnedActivityDock::refresh);
    updateControls();
}

bool LearnedActivityDock::load(const QString& path, QString* error) {
    ASSERT_THREAD(this);
    try {
        require(isEnabled(), QStringLiteral("another Reader operation is running"));
        require(scene_ && protein_ && raw_ && display_, QStringLiteral("open a molecule first"));
        QFile file(path);
        require(file.open(QIODevice::ReadOnly), file.errorString());
        QJsonParseError parseError;
        const auto document = QJsonDocument::fromJson(file.readAll(), &parseError);
        require(parseError.error == QJsonParseError::NoError && document.isObject(), parseError.errorString());
        const auto root = document.object();
        require(indexValue(root.value("schema_version"), "schema_version") == 1, QStringLiteral("expected capture schema 1"));
        require(root.value("coordinate_frame").toString() == QStringLiteral("extraction"),
                QStringLiteral("capture tensors and positions must be in the extraction frame"));
        require(static_cast<std::size_t>(indexValue(root.value("atom_count"), "atom_count")) == protein_->atomCount(),
                QStringLiteral("capture atom count differs from the loaded molecule"));
        const auto name = root.value("model").toString().trimmed();
        require(!name.isEmpty(), QStringLiteral("model name is missing"));
        std::vector<std::size_t> atoms;
        std::set<int> uniqueAtoms;
        for (const auto& value : root.value("atom_indices").toArray()) {
            const int atom = indexValue(value, "atom index");
            require(static_cast<std::size_t>(atom) < protein_->atomCount() && uniqueAtoms.insert(atom).second,
                    QStringLiteral("atom indices must be distinct and in range"));
            atoms.push_back(static_cast<std::size_t>(atom));
        }
        require(!atoms.empty(), QStringLiteral("capture has no atoms"));
        std::map<std::size_t, Frame> frames;
        QStringList names;
        double maximum = 0.0;
        for (const auto& value : root.value("frames").toArray()) {
            const auto object = value.toObject();
            const int original = indexValue(object.value("frame_index"), "frame_index");
            const auto row = raw_->frameRowForOriginalIndex(static_cast<std::size_t>(original));
            require(row.has_value(), QStringLiteral("capture frame %1 is absent from the loaded molecule").arg(original));
            require(frames.count(*row) == 0, QStringLiteral("duplicate capture frame"));
            const auto positions = object.value("positions").toArray();
            require(positions.size() == static_cast<qsizetype>(atoms.size()), QStringLiteral("position count differs from atom indices"));
            for (std::size_t i = 0; i < atoms.size(); ++i) {
                const auto p = numbers(positions[static_cast<qsizetype>(i)], 3);
                const model::Vec3 actual = raw_->atomPosition(*row, atoms[i]);
                require(actual.allFinite() && (actual - model::Vec3(p[0], p[1], p[2])).norm() <= kPositionToleranceA,
                        QStringLiteral("position mismatch at frame %1 atom %2").arg(original).arg(atoms[i]));
            }
            Frame frame;
            QStringList frameNames;
            for (const auto& channelValue : object.value("channels").toArray()) {
                const auto channel = channelValue.toObject();
                const QString channelName = channel.value("name").toString().trimmed();
                require(!channelName.isEmpty() && !frameNames.contains(channelName), QStringLiteral("channel names must be nonempty and distinct"));
                frameNames.push_back(channelName);
                const auto values = channel.value("tensors").toArray();
                require(values.size() == static_cast<qsizetype>(atoms.size()), QStringLiteral("tensor count differs from atom indices"));
                std::vector<model::Mat3> tensors;
                for (const auto& tensorValue : values) {
                    const auto c = numbers(tensorValue, 6);
                    model::Mat3 tensor;
                    tensor << c[0], c[1], c[2], c[1], c[3], c[4], c[2], c[4], c[5];
                    require(std::abs(tensor.trace()) <= kTraceTolerance * std::max(1.0, tensor.norm()),
                            QStringLiteral("hidden 2e tensor is not traceless"));
                    const Eigen::SelfAdjointEigenSolver<model::Mat3> eigen(tensor);
                    require(eigen.info() == Eigen::Success && eigen.eigenvalues().allFinite(), QStringLiteral("tensor eigendecomposition failed"));
                    maximum = std::max(maximum, eigen.eigenvalues().cwiseAbs().maxCoeff());
                    tensors.push_back(tensor);
                }
                frame.channels.push_back(std::move(tensors));
            }
            require(!frameNames.empty(), QStringLiteral("capture has no channels"));
            if (frames.empty()) names = frameNames;
            require(frameNames == names, QStringLiteral("channel names or order change between frames"));
            frames.emplace(*row, std::move(frame));
        }
        require(!frames.empty(), QStringLiteral("capture has no frames"));
        // Commit only after validating the entire capture; a failed load leaves
        // the previous display intact, just like loading a molecule in Reader.
        glyphs_.clear();
        atoms_ = std::move(atoms);
        frames_ = std::move(frames);
        channelNames_ = names;
        modelName_ = name;
        referenceMagnitude_ = maximum > 0.0 ? maximum : 1.0;
        {
            const QSignalBlocker blocker(visible_);
            visible_->setChecked(true);
        }
        for (std::size_t i = 0; i < atoms_.size(); ++i)
            glyphs_.push_back(std::make_unique<TensorGlyphActor>(vtkSmartPointer<vtkRenderer>(scene_->Renderer())));
        updateControls();
        refresh();
        show();
        raise();
        qCInfo(cActivity).noquote() << "capture loaded" << path << "| model=" << modelName_
                                  << "| atoms=" << atoms_.size() << "| frames=" << frames_.size()
                                  << "| channels=" << names.size() << "| reference=" << referenceMagnitude_;
        return true;
    } catch (const std::exception& exception) {
        *error = QString::fromUtf8(exception.what());
        qCWarning(cActivity).noquote() << "capture rejected" << path << "|" << *error;
        return false;
    }
}

bool LearnedActivityDock::configure(const QJsonObject& options, QString* error) {
    ASSERT_THREAD(this);
    try {
        require(!frames_.empty(), QStringLiteral("no activation capture loaded"));
        for (auto it = options.begin(); it != options.end(); ++it)
            require(it.key() == "channel" || it.key() == "radius" || it.key() == "opacity" || it.key() == "visible",
                    QStringLiteral("unknown activity option: ") + it.key());
        int channel = channel_->currentIndex();
        if (options.contains("channel")) channel = indexValue(options.value("channel"), "channel");
        require(channel < channelNames_.size(), QStringLiteral("channel is out of range"));
        const auto setting = [&](const char* name, QDoubleSpinBox* control) {
            if (!options.contains(name)) return control->value();
            const auto value = options.value(name);
            require(value.isDouble() && std::isfinite(value.toDouble()) && value.toDouble() >= control->minimum()
                        && value.toDouble() <= control->maximum(), QString::fromLatin1(name) + " is out of range");
            return value.toDouble();
        };
        const double radius = setting("radius", radius_);
        const double opacity = setting("opacity", opacity_);
        require(!options.contains("visible") || options.value("visible").isBool(), QStringLiteral("visible must be boolean"));
        const QSignalBlocker b1(channel_), b2(radius_), b3(opacity_), b4(visible_);
        channel_->setCurrentIndex(channel);
        radius_->setValue(radius);
        opacity_->setValue(opacity);
        if (options.contains("visible")) visible_->setChecked(options.value("visible").toBool());
        refresh();
        qCInfo(cActivity) << "display configured | channel=" << channel
                         << "radius=" << radius << "opacity=" << opacity
                         << "visible=" << visible_->isChecked();
        return true;
    } catch (const std::exception& exception) {
        *error = QString::fromUtf8(exception.what());
        qCWarning(cActivity).noquote() << "display option rejected |" << *error;
        return false;
    }
}

void LearnedActivityDock::clear() {
    ASSERT_THREAD(this);
    glyphs_.clear();
    atoms_.clear();
    frames_.clear();
    channelNames_.clear();
    modelName_.clear();
    referenceMagnitude_ = 1.0;
    drawnSamples_ = {};
    updateControls();
    if (scene_) scene_->requestRender(MoleculeScene::RenderSource::Overlay);
    qCInfo(cActivity) << "capture cleared";
}

void LearnedActivityDock::updateControls() {
    const QSignalBlocker blocker(channel_);
    channel_->clear();
    channel_->addItems(channelNames_);
    channel_->setEnabled(!frames_.empty());
    radius_->setEnabled(!frames_.empty());
    opacity_->setEnabled(!frames_.empty());
    visible_->setEnabled(!frames_.empty());
    status_->setText(QStringLiteral("No activation capture"));
}

void LearnedActivityDock::refresh() {
    ASSERT_THREAD(this);
    drawnSamples_ = {};
    if (!scene_ || !display_ || !playback_ || frames_.empty()) return;
    const auto frame = static_cast<std::size_t>(playback_->currentFrame());
    const auto found = frames_.find(frame);
    for (auto& glyph : glyphs_) glyph->clear();
    if (found != frames_.end() && visible_->isChecked()) {
        const auto rotation = display_->displayRotation(frame);
        const auto& tensors = found->second.channels[static_cast<std::size_t>(channel_->currentIndex())];
        TensorGlyphActor::Style style;
        style.showSurface = true;
        style.showArrows = false;
        style.surfaceOpacity = opacity_->value();
        style.surfaceReferenceMagnitude = referenceMagnitude_;
        style.ovaloidScale = radius_->value() / 1.5;
        for (std::size_t i = 0; i < atoms_.size(); ++i) {
            const model::Mat3 tensor = rotation * tensors[i] * rotation.transpose();
            const Eigen::SelfAdjointEigenSolver<model::Mat3> eigen(tensor);
            const model::Vec3 center = display_->atomPosition(frame, atoms_[i]);
            const auto& values = eigen.eigenvalues();
            if (values.cwiseAbs().maxCoeff() > 0.0)
                glyphs_[i]->show(center, {values[0], values[1], values[2]}, eigen.eigenvectors(), 0.0, 1.0, style);
            drawnSamples_.append(QJsonObject{{"atom", static_cast<qint64>(atoms_[i])}, {"center", vectorJson(center)},
                {"tensor", tensorJson(tensor)}, {"norm", tensor.norm()},
                {"peak_radius", radius_->value() * values.cwiseAbs().maxCoeff() / referenceMagnitude_}});
        }
    }
    status_->setText(QStringLiteral("%1\n%2 atoms; %3 captured frames\n%4")
        .arg(modelName_).arg(atoms_.size()).arg(frames_.size())
        .arg(found == frames_.end() ? QStringLiteral("No capture for this frame") : channel_->currentText()));
    scene_->requestRender(MoleculeScene::RenderSource::Overlay);
    qCDebug(cActivity) << "frame=" << frame << "channel=" << channel_->currentIndex()
                      << "samples=" << drawnSamples_.size();
}

QJsonObject LearnedActivityDock::state() const {
    ASSERT_THREAD(this);
    return {{"loaded", !frames_.empty()}, {"editable", isEnabled()},
            {"model", modelName_}, {"channels", QJsonArray::fromStringList(channelNames_)},
            {"channel", channel_->currentIndex()}, {"visible", visible_->isChecked()},
            {"radius", radius_->value()}, {"opacity", opacity_->value()}, {"reference_magnitude", referenceMagnitude_},
            {"frame_count", static_cast<qint64>(frames_.size())}, {"samples", drawnSamples_}};
}

void LearnedActivityDock::openCapture() {
    ASSERT_THREAD(this);
    const auto path = QFileDialog::getOpenFileName(this, QStringLiteral("Open activation capture"), {}, QStringLiteral("Activation capture (*.json)"));
    if (path.isEmpty()) return;
    QString error;
    if (!load(path, &error)) QMessageBox::warning(this, QStringLiteral("Activation capture"), error);
}
}  // namespace h5reader::app
