H5 Reader for the external trajectory drive
==========================================

Open "Start Reader.desktop" in 01_Reader/Linux/. Some desktops require one
"Allow launching" action. Reader lists all 176 trajectories and opens their
existing files directly from the drive. It does not copy trajectory datasets.
"Start Reader - Singularity.desktop" selects the installed singularity command.
"Start Reader - Apptainer.desktop" selects installed Apptainer, or offers setup
of the carried Apptainer engine when included.
"Start Reader - bundled Apptainer.desktop" explicitly selects the carried engine,
even when an installed Apptainer is present but does not work.
"Start Reader - compatibility mode.desktop" selects the bundled PRoot runtime.
TRY_NEXT.txt explains these choices. The plain Start Reader.desktop shortcut
automatically chooses installed Apptainer, then installed Singularity, then the
carried engine when included, then the explained PRoot fallback.
ON_THE_DESK.txt gives a short offline guide for the single visit.

Apptainer / Singularity (preferred)
----------------------------------
When installed Apptainer or Singularity is available, the launcher runs the compressed
reader.sif image directly from the external drive. It does not unpack a large
application onto the computer. The image contains its own Ubuntu userspace,
Qt, VTK, software OpenGL and CPU inference libraries; the advisor does not
need to determine the host glibc version or install Qt. The data drive is
mounted read-only inside the container. Only settings, cache, temporary work,
logs and requested output are written to the local workspace.

The carried Apptainer engine is an additional option when neither site engine
is available. Its first-use dialog explains a one-time copy of about 163 MB of
container-engine software to the local workspace. The SIF image and trajectories
stay on the drive. This avoids executing engine binaries from a noexec or
space-containing drive path. The local workspace must contain no whitespace
and allow execution; the carried engine needs glibc 2.28+ and permitted user
namespaces. It was accepted on EL8; an EL7 computer should use site Singularity
or the PRoot alternative. A blocked engine reports its error without silently
installing another runtime. No administrator access or system installation is used.

The default workspace is ~/.local/share/h5reader-portable. To use another SSD:

  H5READER_WORKSPACE=/absolute/path/on/ssd /bin/sh start-reader.sh

Workspace and source drive must be separate, non-overlapping directories.
The image uses CPU rendering and CPU inference. GPU devices are not requested.
Requirements: x86_64 Linux, X11 or XWayland, a functioning Apptainer 1.0+ or
Singularity 3.7+ installation, and site policy permitting the required binds.
The launcher explicitly disables NVIDIA/ROCm integration, including site defaults.
CentOS 7 with Singularity 3.8.4 and AlmaLinux 8 with Singularity 3.8.4 and the
carried Apptainer have passed acceptance; see ACCEPTANCE.txt. The advisor's
actual kernel, site policy and desktop still need the compatibility check.

Checking the advisor's computer
-------------------------------
Open "Check compatibility.desktop", or run:

  /bin/sh /path/to/01_Reader/Linux/check-compatibility.sh

Add --backend singularity or --backend apptainer to check that exact command.

This quick check runs the actual image and Reader's version-only Qt/X11 startup.
It can also use a previously prepared carried engine; it never copies that
engine during the check. If local engine setup is needed it explains which
launcher to use first.
It uses tiny synthetic files to verify a read-only source bind and a writable
local workspace. It does not open trajectories, unpack the fallback application,
install anything or send the report anywhere. Its report is saved under the local
workspace's state/ directory. --workspace /absolute/local/path selects another
local workspace; --report /absolute/local/file.txt chooses the report filename.
The report records host/kernel, distribution, host glibc, container version,
display availability, each result and the exact runtime/Qt error if one occurs.

Host glibc does not need to match the application's glibc: the image contains
its own. The host kernel, container installation and site policy still matter.
The check exercises these together instead of guessing from a distribution name.
It does not replace a full trajectory/rendering/inference acceptance session.
If a site provides Singularity through environment modules, load its module in
a terminal and run the check from that terminal. Profiles and modules are never
loaded automatically. A terminal-only session without DISPLAY can check the
container and mounts; Qt/X11 must then be checked from the desktop session.

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

Open "Start Reader - compatibility mode.desktop" or use --backend proot to
request the fallback. Use --backend apptainer to require Apptainer.
Container launch failures are reported; they do not silently trigger
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

To include an independently tested relocatable Apptainer tree, add
--bundled-apptainer /staging/tested-apptainer-tree. It is carried separately from
the image; its file hashes, symlinks and modes are recorded in runtime/.

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
  https://apptainer.org/docs/admin/latest/installation.html
  https://doc.qt.io/qt-6.10/supported-platforms.html
  https://proot-me.github.io/
  https://specifications.freedesktop.org/desktop-entry/latest-single/
