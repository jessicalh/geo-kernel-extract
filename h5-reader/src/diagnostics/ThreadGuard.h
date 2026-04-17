// ThreadGuard — assert thread affinity at method entry.
//
// ASSERT_THREAD(obj) — at the top of any method with thread-affinity
// expectations. If the current thread differs from obj->thread(), log a
// critical message and (in Debug builds) abort via Q_ASSERT. Use in
// slots that must run on the GUI thread, in worker methods that must not
// touch GUI state, anywhere a thread violation would cause undefined
// behaviour.

#pragma once

#include <QLoggingCategory>
#include <QObject>
#include <QThread>

namespace h5reader::diagnostics {

Q_DECLARE_LOGGING_CATEGORY(cThreadGuard)

}  // namespace

#define ASSERT_THREAD(obj) \
    do { \
        const QObject* _tg_obj = (obj); \
        if (_tg_obj && QThread::currentThread() != _tg_obj->thread()) { \
            qCCritical(::h5reader::diagnostics::cThreadGuard).nospace() \
                << "[thread] affinity violation at " \
                << __FILE__ << ":" << __LINE__ \
                << " — current=" << QThread::currentThread() \
                << " expected=" << _tg_obj->thread() \
                << " object=" << _tg_obj->metaObject()->className(); \
            Q_ASSERT_X(false, "ASSERT_THREAD", "Thread affinity violation"); \
        } \
    } while (0)
