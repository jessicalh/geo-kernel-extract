#pragma once

#include <QtGlobal>

#ifdef Q_OS_MACOS
#include <QAbstractItemView>
#include <QPainter>
#include <QProxyStyle>
#include <QStyleOptionToolButton>
#include <QToolBar>
#include <QToolButton>

namespace h5reader::app {

// Local corrections to Qt 6.12's native Mac style. Recheck both drawing
// paths when changing Qt versions; leave the application style untouched.
class MacWidgetStyle final : public QProxyStyle {
public:
    explicit MacWidgetStyle(const QString& styleName)
        : QProxyStyle(styleName) {}

    void drawPrimitive(PrimitiveElement element, const QStyleOption* option,
                       QPainter* painter, const QWidget* widget = nullptr) const override {
        if (element == PE_IndicatorItemViewItemCheck) {
            // QMacStyle passes the item-view rectangle through to AppKit as
            // both a native view frame and its drawing bounds. With a nonzero
            // origin, the checkbox disappears on macOS 27 / Qt 6.12. Draw the
            // same native indicator in local coordinates; keep its original
            // size, clipping, state, palette, and the delegate's hit rectangle.
            QStyleOption local(*option);
            painter->save();
            painter->translate(local.rect.topLeft());
            local.rect.moveTopLeft(QPoint(0, 0));
            QProxyStyle::drawPrimitive(element, &local, painter, widget);
            painter->restore();
            return;
        }
        QProxyStyle::drawPrimitive(element, option, painter, widget);
    }

    void drawControl(ControlElement element, const QStyleOption* option,
                     QPainter* painter, const QWidget* widget = nullptr) const override {
        if (element == CE_ToolButtonLabel) {
            // The native toolbar background is translucent, but QMacStyle
            // uses selected-menu text for checked text buttons. The common
            // label renderer respects the button palette and disabled state.
            const auto* button = qstyleoption_cast<const QStyleOptionToolButton*>(option);
            if (button && !(button->features & QStyleOptionToolButton::Arrow)
                && (button->toolButtonStyle == Qt::ToolButtonTextOnly
                    || button->icon.isNull())) {
                QCommonStyle::drawControl(element, option, painter, widget);
                return;
            }
        }
        QProxyStyle::drawControl(element, option, painter, widget);
    }
};

// Call after adding the toolbar's actions/widgets. Use a separately owned
// native style: wrapping QApplication::style() would take its ownership and
// affect unrelated controls. No palette, stylesheet, or SDK changes.
inline void configureMacToolbarText(QToolBar* toolbar) {
    const QString styleName = toolbar->style()->name();
    if (styleName.compare(QStringLiteral("macos"), Qt::CaseInsensitive) != 0)
        return;

    auto* style = new MacWidgetStyle(styleName);
    style->setParent(toolbar);
    for (auto* button : toolbar->findChildren<QToolButton*>(
             QString(), Qt::FindDirectChildrenOnly)) {
        if (button->icon().isNull() && !button->text().isEmpty())
            button->setStyle(style);
    }
}

inline void configureMacItemViewChecks(QAbstractItemView* view) {
    const QString styleName = view->style()->name();
    if (styleName.compare(QStringLiteral("macos"), Qt::CaseInsensitive) != 0)
        return;

    auto* style = new MacWidgetStyle(styleName);
    style->setParent(view);
    view->setStyle(style);
}

} // namespace h5reader::app
#endif // Q_OS_MACOS
