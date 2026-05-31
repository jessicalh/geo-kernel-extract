// ChargeStore — FF14SB topol.top charge reader for rediscover multipoles.

#pragma once

#include <QString>

namespace h5reader::model {
class QtProtein;
}

namespace h5reader::rediscover {

bool LoadFf14sbChargesFromTopol(const QString& topolPath, model::QtProtein& protein,
                                QString* err_out = nullptr);

}  // namespace h5reader::rediscover
