#ifndef PROGRESS_REPORT_H
#define PROGRESS_REPORT_H

/** Optional hook so a long solver can report done/total while it is still running.
 *  The default is a no-op. The WASM page installs a hook that posts to the UI. */
using ProgressHook = void (*)(int done, int total);

inline ProgressHook& progressHook() {
    static ProgressHook hook = nullptr;
    return hook;
}

inline int& progressBucket() {
    static int bucket = -1;
    return bucket;
}

inline void progressReset() {
    progressBucket() = -1;
}

inline void progressReport(int done, int total) {
    if (!progressHook() || total < 1) return;
    if (done < 0) done = 0;
    if (done > total) done = total;
    const int buckets = total < 40 ? total : 100;
    const int bucket = (done * buckets) / total;
    if (bucket == progressBucket() && done != total) return;
    progressBucket() = bucket;
    progressHook()(done, total);
}

#endif
