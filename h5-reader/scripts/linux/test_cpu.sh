#!/usr/bin/env bash
# Run configured native and REST acceptance tests using Mesa software rendering.
set -euo pipefail
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 BUILD_DIRECTORY [ctest arguments...]" >&2
    exit 2
fi
build_directory=$1
shift
exec env -u LD_LIBRARY_PATH -u LD_PRELOAD \
    -u QT_PLUGIN_PATH -u QT_QPA_PLATFORM_PLUGIN_PATH \
    CUDA_VISIBLE_DEVICES=-1 HIP_VISIBLE_DEVICES=-1 ROCR_VISIBLE_DEVICES=-1 \
    H5READER_EXPERIMENTAL_SHIELDING_ML_DEVICE=cpu \
    LIBGL_ALWAYS_SOFTWARE=1 GALLIUM_DRIVER=llvmpipe \
    __GLX_VENDOR_LIBRARY_NAME=mesa QT_QPA_PLATFORM=xcb \
    QT_OPENGL=software QT_QUICK_BACKEND=software \
    QT_FFMPEG_ENCODING_HW_DEVICE_TYPES=, QT_FFMPEG_DECODING_HW_DEVICE_TYPES=, \
    QT_DISABLE_HW_TEXTURES_CONVERSION=1 OMP_NUM_THREADS=2 \
    xvfb-run -a -s '-screen 0 1280x720x24 +extension GLX' \
    ctest --test-dir "$build_directory" --output-on-failure "$@"
