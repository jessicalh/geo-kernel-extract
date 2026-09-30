H5 Reader for the external trajectory drive
==========================================

Open "Start Reader.desktop" in 01_Reader/Linux/. Some desktops require one
"Allow launching" action. Reader lists all 176 trajectories and opens their
existing files directly from the drive. It does not copy trajectory datasets.

Apptainer / Singularity (preferred)
----------------------------------
When Apptainer or Singularity is available, the launcher runs the compressed
reader.sif image directly from the external drive. It does not unpack a large
application onto the computer. The image contains its own Ubuntu userspace,
Qt, VTK, software OpenGL and CPU inference libraries; the advisor does not
need to determine the host glibc version or install Qt. The data drive is
mounted read-only inside the container. Only settings, cache, temporary work,
logs and requested output are written to the local workspace.

The default workspace is ~/.local/share/h5reader-portable. To use another SSD:

  H5READER_WORKSPACE=/absolute/path/on/ssd /bin/sh start-reader.sh

Workspace and source drive must be separate, non-overlapping directories.
The image uses CPU rendering and CPU inference. GPU devices are not requested.
Requirements: x86_64 Linux, X11 or XWayland, a functioning Apptainer/Singularity
installation, and site policy permitting the required binds. The host kernel
and container setup still require acceptance on the advisor's actual machine;
no particular RHEL or CentOS installation is claimed tested in advance.

Optional PRoot fallback
-----------------------
Reader is on the external drive. If Apptainer/Singularity is absent, the
computer needs a compatible application runtime on local storage. A first-time
dialog explains this optional one-time setup: unpack about 1.3 GB of application
software locally. Choose "Set up Reader" or "Cancel". A progress dialog can
cancel extraction. Later launches reuse the prepared application. Trajectories
stay on the drive and open directly. This fallback needs ptrace permission
and local storage that permits program execution. PRoot is a compatibility
runtime, not a security sandbox, and cannot enforce a read-only source mount.

Use --backend proot to request the fallback, or --backend apptainer to require
Apptainer. Container launch failures are reported; they do not silently trigger
a large fallback installation. Fallback setup needs zenity, kdialog, xmessage,
or an interactive terminal to obtain confirmation. --non-interactive-setup is
an explicit developer acceptance option only.

The shell launcher can be interpreted with /bin/sh from a noexec source mount.
The SIF is read as data; fallback executables run from local storage. Older
"gio launch" versions have an upstream bug that loses the desktop file's %k
location. Normal activation through Gio.DesktopAppInfo.new_from_filename is
tested separately. On such a command line use /bin/sh /path/to/start-reader.sh.

Packaging on the development Ubuntu 24.04 host
--------------------------------------------
Download these inputs into writable SSD staging:

  https://cdimage.ubuntu.com/ubuntu-base/releases/24.04/release/ubuntu-base-24.04.4-base-amd64.tar.gz
  https://proot.gitlab.io/proot/bin/proot

build_runtime.py checks the pinned Ubuntu SHA256 and records the PRoot version
and SHA256, Reader binary hash, local catalog hash/count, SIF hash and fallback
archive hash. The installed bundle must contain linux-local-library.json beside
h5reader. Software Mesa dependencies are copied from the packaging host; no
host GPU computation or global install is needed.

  python3 build_runtime.py --bundle /path/to/installed-reader \
    --base /staging/ubuntu-base-24.04.4-base-amd64.tar.gz \
    --proot /staging/proot --output /staging/mock-drive/01_Reader/Linux \
    --build-root /staging/temporary-rootfs

Apptainer converts the prepared userspace to a SIF offline. The deliverable
contains the launcher, image, optional fallback archive and manifest. Build
directories stay separate. Neither packaging nor acceptance installs onto the
actual Provenance drive; that remains an explicitly coordinated later step.

Acceptance uses Reader's built-in REST and snapshot tools:

  python3 smoke_runtime.py --launcher /staging/mock-drive/01_Reader/Linux/start-reader.sh \
    --source-root /staging/mock-drive --workspace /staging/acceptance-workspace \
    --fixture 01_Reader/datasets/SinglePose/example/run.LGS --check-ml \
    --check-local-library-key local-bmr7057 --output /staging/acceptance

--check-ml requires complete F006 inputs. The test exercises rendered snapshots,
CPU prediction, all 176 catalog entries, direct opening and frame navigation,
no trajectory copies, graceful shutdown, and scientific-file SHA256 equality
before and after the session. Xvfb is used if no display is supplied. Additional
Reader arguments may be passed after the shell launcher's -- separator.

References:
  https://apptainer.org/docs/user/latest/build_a_container.html
  https://apptainer.org/docs/user/latest/bind_paths_and_mounts.html
  https://proot-me.github.io/
  https://specifications.freedesktop.org/desktop-entry/latest-single/
