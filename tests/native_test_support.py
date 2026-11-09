"""Native test discovery without silently accepting a broken explicit override."""

import os
from pathlib import Path
import shutil
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REQUIRE_NATIVE_ENV = "SJARACNE_REQUIRE_NATIVE_TESTS"
_MACHO_MAGICS = {
    bytes.fromhex(value)
    for value in (
        "feedface", "cefaedfe", "feedfacf", "cffaedfe",
        "cafebabe", "bebafeca", "cafebabf", "bfbafeca",
    )
}


def _runtime_problem(path):
    """Check permissions and file format, without invoking the executable."""
    if not path.is_file():
        return "not a file"
    try:
        with path.open("rb") as handle:
            magic = handle.read(4)
    except OSError as error:
        return "cannot read file: {}".format(error)

    if sys.platform == "win32":
        if not magic.startswith(b"MZ"):
            return "not a Windows executable (missing MZ header)"
    else:
        if not os.access(str(path), os.X_OK):
            return "not executable (check file permissions)"
        # A shebang permits explicit executable wrappers as well as native files.
        if sys.platform.startswith("linux") and not (
                magic == b"\x7fELF" or magic.startswith(b"#!")):
            return "not a Linux executable (expected ELF or a shebang wrapper)"
        if sys.platform == "darwin" and not (
                magic in _MACHO_MAGICS or magic.startswith(b"#!")):
            return "not a macOS executable (expected Mach-O or a shebang wrapper)"
    return None


def find_native_executable(name="sjaracne.exe", override_env="SJARACNE_TEST_EXE",
                           project_root=None):
    """Return an absolute executable path, or None when native tests are optional.

    An explicit override is authoritative: an invalid path never falls back to
    a different binary. Required mode makes missing/unusable auto-discovered
    binaries an error instead of allowing each test class to skip.
    """
    if override_env in os.environ:
        configured = os.environ[override_env]
        if not configured:
            raise RuntimeError("{} is set but empty".format(override_env))
        path = Path(configured).resolve()
        problem = _runtime_problem(path)
        if problem:
            raise RuntimeError("{}={} is unusable: {}".format(
                override_env, path, problem))
        return str(path)

    root = PROJECT_ROOT if project_root is None else Path(project_root)
    candidates = [root / "SJARACNe" / "bin" / name, shutil.which(name)]
    problems = []
    for candidate in candidates:
        if candidate is None:
            continue
        path = Path(candidate).resolve()
        problem = _runtime_problem(path)
        if problem is None:
            return str(path)
        problems.append("{} ({})".format(path, problem))

    if os.environ.get(REQUIRE_NATIVE_ENV) == "1":
        raise RuntimeError(
            "{}=1 requires {}; build it or set {}. Checked: {}".format(
                REQUIRE_NATIVE_ENV, name, override_env, "; ".join(problems)))
    return None
