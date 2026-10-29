#!/usr/bin/env python3
"""Build a clean source archive and verify its wheel in an isolated environment.

This is an explicit packaging/network-install check, not part of unit-test
discovery. Prepare the locked regression group (or pip build requirements)
first. Build and installation
logs, archives and fingerprints are retained; temporary source and virtual
environment copies are removed automatically. No source checkout is installed
into the fresh environment and no dataset benchmark jobs are launched.
"""

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import venv
import zipfile


ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ASSETS = (
    'SJARACNe/__init__.py', 'SJARACNe/sjaracne.py',
    'SJARACNe/apmi_cutoff.py', 'SJARACNe/workflow_options.py',
    'SJARACNe/bin/sjaracne.exe', 'SJARACNe/bin/apmi_null_generator.exe',
    'SJARACNe/bin/QC_input.py', 'SJARACNe/bin/create_consensus_network.py',
    'SJARACNe/bin/wdl_workflow.py', 'SJARACNe/config/apmi_pchip.json',
    'SJARACNe/config/config_threshold.txt', 'SJARACNe/config/config_kernel.txt',
    'SJARACNe/config/TF_list.txt', 'SJARACNe/config/SIG_list.txt',
    'SJARACNe/cwl/sjaracne.cwl', 'SJARACNe/cwl/sjaracne_workflow.cwl',
    'SJARACNe/wdl/sjaracne.wdl', 'SJARACNe/wdl/sprocket.local.toml',
)
SOURCE_ASSETS = (
    'setup.py', 'pyproject.toml', 'uv.lock', '.python-version', 'MANIFEST.in', 'version.py',
    'SJARACNe/Makefile', 'SJARACNe/src/main.cpp', 'SJARACNe/src/apmi.cpp',
    'SJARACNe/src/apmi_pchip_data.h', 'tools/apmi_null_generator.cpp',
    'docs/BACKWARD_COMPATIBILITY.md',
) + tuple(name for name in PACKAGE_ASSETS if not name.endswith('.exe'))


PROBE = r'''
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import SJARACNe

package = Path(SJARACNe.__file__).resolve().parent
prefix = Path(sys.prefix).resolve()
assert os.path.commonpath((str(package), str(prefix))) == str(prefix), str(package)
assets = json.loads(sys.argv[2])
hashes = {}
for name in assets:
    path = package.parent / name
    assert path.is_file(), name
    hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
native = package / 'bin/sjaracne.exe'
assert os.access(native, os.X_OK), str(native)
native_help = subprocess.run([str(native), '--help'],
                             capture_output=True, text=True)
assert native_help.returncode in (0, 1), native_help.stderr
assert 'Usage:' in native_help.stdout
command = Path(sys.executable).parent / 'sjaracne'
cli_help = subprocess.run([str(command), 'local', '--help'], check=True,
                          capture_output=True, text=True)
assert '--network-type' in cli_help.stdout
cutoff = SJARACNe.get_mi_cutoff(80, 0.001)
assert math.isfinite(cutoff) and cutoff > 0, cutoff
record = {
    'python': sys.version, 'prefix': str(prefix), 'package_root': str(package),
    'versions': {name: importlib.metadata.version(name)
                 for name in ('SJARACNe', 'numpy', 'pandas', 'scipy', 'cwltool')},
    'assets_sha256': hashes, 'pchip_m80_p001': cutoff,
    'native_help': native_help.stdout, 'native_help_returncode': native_help.returncode,
    'cli_help': cli_help.stdout,
}
Path(sys.argv[1]).write_text(json.dumps(record, indent=2, sort_keys=True) + '\n',
                           encoding='utf-8')
'''


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n',
                          encoding='utf-8')


def check_members(names, expected):
    missing = sorted(set(expected) - set(names))
    if missing:
        raise RuntimeError('Package is missing required assets: ' + ', '.join(missing))
    retired = sorted(name for name in names
                     if 'cwlexec' in name.lower() or name.startswith('examples/lsf/'))
    if retired:
        raise RuntimeError('Retired runner assets remain packaged: ' + ', '.join(retired))


def inspect_source_archive(path):
    with tarfile.open(path, 'r:gz') as archive:
        names = [member.name.split('/', 1)[1] for member in archive.getmembers()
                 if member.isfile() and '/' in member.name]
    check_members(names, SOURCE_ASSETS)
    if any(name.endswith('.exe') for name in names):
        raise RuntimeError('Clean source archive unexpectedly contains prebuilt executables')
    if 'BACKWARD_COMPATIBILITY.md' in names:
        raise RuntimeError('Old root compatibility guide remains packaged')


def inspect_wheel(path):
    with zipfile.ZipFile(path) as wheel:
        names = wheel.namelist()
        check_members(names, PACKAGE_ASSETS)
        metadata = [name for name in names if name.endswith('.dist-info/WHEEL')]
        if len(metadata) != 1:
            raise RuntimeError('Expected one wheel metadata record')
        record = wheel.read(metadata[0]).decode('utf-8')
    if 'Root-Is-Purelib: false' not in record:
        raise RuntimeError('Native wheel incorrectly advertises pure Python contents')
    tags = [line.split(': ', 1)[1] for line in record.splitlines() if line.startswith('Tag: ')]
    if not tags or any(tag.endswith('-any') for tag in tags):
        raise RuntimeError('Native wheel incorrectly advertises platform independence')
    return tags


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--requirements', type=Path,
                        default=ROOT / 'tests/requirements-ci.txt',
                        help='Runtime constraints, preferably exported from uv.lock')
    parser.add_argument('--build-requirements', type=Path,
                        default=ROOT / 'tests/requirements-build.txt',
                        help='Build-tool requirements; may use the same locked export')
    arguments = parser.parse_args()
    output = arguments.output_dir.expanduser().resolve()
    requirements = arguments.requirements.expanduser().resolve()
    build_requirements = arguments.build_requirements.expanduser().resolve()
    if not requirements.is_file():
        parser.error('Runtime requirements file does not exist: ' + str(requirements))
    if not build_requirements.is_file():
        parser.error('Build requirements file does not exist: ' + str(build_requirements))
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error('Output must be a new or empty directory: ' + str(output))
    output.mkdir(parents=True, exist_ok=True)
    record = {
        'schema': 'sjaracne-source-install-v1', 'complete': False,
        'started_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'python': sys.version, 'requirements': str(requirements),
        'requirements_sha256': sha256(requirements),
        'build_requirements': str(build_requirements),
        'build_requirements_sha256': sha256(build_requirements), 'steps': [],
    }
    manifest = output / 'manifest.json'
    write_json(manifest, record)
    environment = os.environ.copy()
    environment.pop('PYTHONPATH', None)
    environment.pop('PYTHONHOME', None)
    environment.update(PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1',
                       PIP_DISABLE_PIP_VERSION_CHECK='1')

    def run(name, command, cwd):
        print('Source-install check: ' + name, flush=True)
        log = output / (name + '.log')
        with log.open('w', encoding='utf-8') as handle:
            result = subprocess.run([str(part) for part in command], cwd=cwd,
                                    env=environment, stdout=handle,
                                    stderr=subprocess.STDOUT, text=True)
        record['steps'].append({'name': name, 'argv': [str(part) for part in command],
                                'returncode': result.returncode,
                                'log': log.name, 'log_sha256': sha256(log)})
        write_json(manifest, record)
        if result.returncode:
            raise RuntimeError(name + ' failed; see ' + str(log))
        return log

    try:
        # Build in a disposable copy, with no cached native binaries or objects.
        # This is safe to run concurrently for different Python versions.
        with tempfile.TemporaryDirectory(prefix='sjaracne-source-install-') as temporary:
            scratch = Path(temporary)
            source = scratch / 'source'
            shutil.copytree(ROOT, source, ignore=shutil.ignore_patterns(
                '.git', '.venv', 'venv', '__pycache__', '*.pyc', '*.egg-info',
                'build', 'dist', '.pytest_cache', 'results', '*.exe'))
            run('build-sdist', [sys.executable, '-m', 'build', '--sdist',
                               '--no-isolation', '--outdir', output / 'dist'], source)
            archives = sorted((output / 'dist').glob('*.tar.gz'))
            if len(archives) != 1:
                raise RuntimeError('Expected exactly one source archive')
            archive = archives[0]
            inspect_source_archive(archive)
            run('build-wheel-from-sdist', [sys.executable, '-m', 'pip', 'wheel',
                                         '--no-cache-dir', '--no-deps',
                                         '--no-build-isolation', '--wheel-dir',
                                         output / 'wheels', archive], scratch)
            wheels = sorted((output / 'wheels').glob('*.whl'))
            if len(wheels) != 1:
                raise RuntimeError('Expected exactly one built wheel')
            wheel = wheels[0]
            tags = inspect_wheel(wheel)
            fresh = scratch / 'fresh-env'
            venv.EnvBuilder(with_pip=True).create(fresh)
            python = fresh / ('Scripts/python.exe' if os.name == 'nt' else 'bin/python')
            run('fresh-build-tools', [python, '-m', 'pip', 'install', '-r',
                                      build_requirements], scratch)
            run('fresh-wheel-install', [python, '-m', 'pip', 'install', '--no-cache-dir',
                                        '--constraint', requirements,
                                        str(wheel) + '[cwl]'], scratch)
            run('fresh-pip-check', [python, '-m', 'pip', 'check'], scratch)
            run('fresh-pip-freeze', [python, '-m', 'pip', 'freeze', '--all'], scratch)
            run('installed-runtime', [python, '-c', PROBE,
                                      output / 'installed-runtime.json',
                                      json.dumps(PACKAGE_ASSETS)], scratch)
            record.update(source_archive={'name': str(archive.relative_to(output)),
                                          'sha256': sha256(archive)},
                          wheel={'name': str(wheel.relative_to(output)),
                                 'sha256': sha256(wheel), 'tags': tags})
        record['complete'] = True
        print('Clean source archive, native wheel and fresh installation passed.', flush=True)
        return 0
    except Exception as error:
        record['error'] = str(error)
        print(str(error), file=sys.stderr)
        return 1
    finally:
        record['finished_at_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        write_json(manifest, record)


if __name__ == '__main__':
    sys.exit(main())
