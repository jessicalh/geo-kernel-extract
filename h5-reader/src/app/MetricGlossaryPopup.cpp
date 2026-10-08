#include "MetricGlossaryPopup.h"

#include "../model/MetricGlossary.h"

#include <QApplication>
#include <QFrame>
#include <QGuiApplication>
#include <QLabel>
#include <QScreen>
#include <QVBoxLayout>

#include <optional>

namespace h5reader::app {

namespace {

QLabel* addHeading(QVBoxLayout* layout, const QString& text, QWidget* parent) {
    auto* label = new QLabel(text, parent);
    label->setTextFormat(Qt::PlainText);
    QFont font = label->font();
    font.setBold(true);
    label->setFont(font);
    layout->addWidget(label);
    return label;
}

QLabel* addBody(QVBoxLayout* layout, const QString& objectName,
                const QString& text, QWidget* parent) {
    auto* label = new QLabel(text, parent);
    label->setTextFormat(Qt::PlainText);
    label->setObjectName(objectName);
    label->setWordWrap(true);
    label->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(label);
    return label;
}

class MetricGlossaryPopup final : public QFrame {
public:
    MetricGlossaryPopup(const QString& titleText,
                        const model::MetricGlossaryEntry& entry,
                        QWidget* parent)
        : QFrame(parent, Qt::Popup) {
        setObjectName(QStringLiteral("MetricGlossaryPopup"));
        setAttribute(Qt::WA_DeleteOnClose);
        setFrameShape(QFrame::StyledPanel);
        setMaximumWidth(520);

        auto* layout = new QVBoxLayout(this);
        layout->setContentsMargins(14, 12, 14, 12);
        layout->setSpacing(5);

        auto* title = new QLabel(titleText, this);
        title->setTextFormat(Qt::PlainText);
        title->setObjectName(QStringLiteral("glossaryTitle"));
        QFont titleFont = title->font();
        titleFont.setBold(true);
        titleFont.setPointSizeF(titleFont.pointSizeF() + 1.0);
        title->setFont(titleFont);
        title->setWordWrap(true);
        title->setTextInteractionFlags(Qt::TextSelectableByMouse);
        layout->addWidget(title);

        layout->addSpacing(3);
        addHeading(layout, QStringLiteral("Meaning"), this);
        addBody(layout, QStringLiteral("glossaryMeaning"), entry.meaning, this);
        layout->addSpacing(3);
        addHeading(layout, QStringLiteral("Calculation"), this);
        addBody(layout, QStringLiteral("glossaryCalculation"), entry.calculation, this);
        layout->addSpacing(3);
        addHeading(layout, QStringLiteral("Origin"), this);
        addBody(layout, QStringLiteral("glossaryOrigin"), entry.origin, this);
    }

    void showAt(const QPoint& globalPosition) {
        adjustSize();
        QPoint position = globalPosition + QPoint(6, 6);
        QScreen* screen = QGuiApplication::screenAt(globalPosition);
        if (!screen)
            screen = QApplication::primaryScreen();
        if (screen) {
            const QRect available = screen->availableGeometry();
            position.setX(qBound(available.left(), position.x(),
                                 available.right() - width() + 1));
            position.setY(qBound(available.top(), position.y(),
                                 available.bottom() - height() + 1));
        }
        move(position);
        show();
        setFocus(Qt::PopupFocusReason);
    }
};

}  // namespace

void ShowMetricGlossaryPopup(const model::SignalDescriptor& descriptor,
                             const QPoint& globalPosition,
                             QWidget* parent) {
    const std::optional<model::MetricGlossaryEntry> entry = model::MetricGlossaryFor(descriptor);
    Q_ASSERT_X(entry.has_value(), "ShowMetricGlossaryPopup",
               "Every catalog descriptor must have a glossary entry");
    if (!entry)
        return;
    ShowMetricGlossaryPopup(descriptor.label, *entry, globalPosition, parent);
}

void ShowMetricGlossaryPopup(const QString& title,
                             const model::MetricGlossaryEntry& entry,
                             const QPoint& globalPosition,
                             QWidget* parent) {
    auto* popup = new MetricGlossaryPopup(title, entry, parent);
    popup->showAt(globalPosition);
}

}  // namespace h5reader::app
