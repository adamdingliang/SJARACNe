#!/usr/bin/env python3
"""Probe automatic editable/wheel native builds, offline reuse and cache refresh.

This explicit Linux/WSL integration check is separate from unit-test discovery.
It uses only a disposable source copy and virtual environments, preserving the
checkout and any existing environment. The supplied cache must be new/empty;
it and all evidence remain available after the temporary copies are removed.
No biological dataset or resampling benchmark is run.
"""

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
RUNTIME_PROBE = r'''
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

import SJARACNe

package = Path(SJARACNe.__file__).resolve().parent
prefix = Path(sys.prefix).resolve()
if sys.argv[2] == 'installed':
    assert os.path.commonpath((str(package), str(prefix))) == str(prefix), str(package)
else:
    assert package == Path(sys.argv[2]).resolve(), str(package)
assert shutil.which('uv') is None, 'Runtime probe must not find uv on PATH'
native = package / 'bin/sjaracne.exe'
native_hashes = {}
for filename in ('sjaracne.exe', 'apmi_null_generator.exe'):
    binary = package / 'bin' / filename
    assert binary.is_file() and os.access(binary, os.X_OK), str(binary)
    native_hashes[filename] = hashlib.sha256(binary.read_bytes()).hexdigest()
help_result = subprocess.run([str(native), '--help'], capture_output=True, text=True)
assert help_result.returncode in (0, 1), help_result.stderr
assert 'Usage:' in help_result.stdout
command = Path(sys.executable).parent / 'sjaracne'
cli = subprocess.run([str(command), 'local', '--help'], check=True,
                     capture_output=True, text=True)
assert '--network-type' in cli.stdout
cutoff = SJARACNe.get_mi_cutoff(80, 0.001)
assert math.isfinite(cutoff) and cutoff > 0
record = {
    'python': sys.version, 'package_root': str(package),
    'versions': {name: importlib.metadata.version(name)
                 for name in ('SJARACNe', 'numpy', 'pandas', 'scipy')},
    'uv_available_on_path': shutil.which('uv'),
    'pchip_m80_p001': cutoff,
    'native_sha256': hashlib.sha256(native.read_bytes()).hexdigest(),
    'native_binaries_sha256': native_hashes,
    'config_asset': 'config/TF_list.txt',
    'config_sha256': hashlib.sha256((package / 'config/TF_list.txt').read_bytes()).hexdigest(),
    'native_help_returncode': help_result.returncode,
}
Path(sys.argv[1]).write_text(json.dumps(record, indent=2, sort_keys=True) + '\n',
                           encoding='utf-8')
'''


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, record):
    Path(path).write_text(json.dumps(record, indent=2, sort_keys=True) + '\n',
                          encoding='utf-8')


def package_was_built(log):
    return bool(re.search(r'\b(?:Building|Built)\s+[`\'"]?sjaracne(?:\s|@|[`\'"])',
                          log, re.IGNORECASE))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--uv', required=True, type=Path, help='Exact uv 0.11.16 executable')
    parser.add_argument('--python', required=True, type=Path,
                        help='Existing Linux Python 3.11 or 3.14 interpreter')
    parser.add_argument('--cache-dir', required=True, type=Path,
                        help='New/empty cache reserved for this probe; retained afterwards')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='New/empty evidence directory')
    args = parser.parse_args()
    if platform.system() != 'Linux':
        parser.error('Use Linux or WSL for native build and uv deployment checks')
    uv = args.uv.expanduser().absolute()
    # Do not resolve interpreter symlinks: a venv launcher may be intentional.
    python = args.python.expanduser().absolute()
    for name, executable in (('uv', uv), ('Python', python)):
        if not executable.is_file() or not os.access(executable, os.X_OK):
            parser.error(name + ' executable is missing: ' + str(executable))
    output = args.output_dir.expanduser().resolve()
    cache = args.cache_dir.expanduser().resolve()
    for name, path in (('Output', output), ('Cache', cache)):
        if path.exists() and (not path.is_dir() or any(path.iterdir())):
            parser.error(name + ' must be a new or empty directory: ' + str(path))
    if output == cache or output in cache.parents or cache in output.parents:
        parser.error('Output and cache directories must be separate, not nested')
    output.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    manifest = output / 'manifest.json'
    record = {
        'schema': 'sjaracne-uv-install-v2', 'complete': False,
        'started_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'uv': str(uv), 'uv_sha256': sha256(uv), 'python': str(python),
        'cache': str(cache), 'source_lock_sha256': sha256(ROOT / 'uv.lock'), 'steps': [],
    }
    write_json(manifest, record)
    environment = os.environ.copy()
    for key in ('PYTHONPATH', 'PYTHONHOME', 'VIRTUAL_ENV', 'UV_PROJECT_ENVIRONMENT', 'UV_PYTHON',
                'UV_OFFLINE', 'UV_NO_BUILD', 'UV_NO_SYNC', 'UV_NO_MANAGED_PYTHON'):
        environment.pop(key, None)
    environment.update(UV_CACHE_DIR=str(cache), UV_PYTHON_DOWNLOADS='never',
                       UV_LINK_MODE='copy', UV_NO_PROGRESS='1',
                       PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1')

    def run(name, command, cwd, env=None, expect_failure=False):
        print('uv install check: ' + name, flush=True)
        path = output / (name + '.log')
        with path.open('w', encoding='utf-8') as handle:
            result = subprocess.run([str(item) for item in command], cwd=cwd,
                                    env=environment if env is None else env,
                                    stdout=handle, stderr=subprocess.STDOUT, text=True)
        content = path.read_text(encoding='utf-8')
        record['steps'].append({'name': name, 'argv': [str(item) for item in command],
                                'returncode': result.returncode, 'log': path.name,
                                'expected_failure': expect_failure,
                                'passed': bool(result.returncode) == expect_failure,
                                'log_sha256': sha256(path),
                                'sjaracne_build_reported': package_was_built(content)})
        write_json(manifest, record)
        if expect_failure and not result.returncode:
            raise RuntimeError(name + ' unexpectedly succeeded; see ' + str(path))
        if result.returncode and not expect_failure:
            raise RuntimeError(name + ' failed; see ' + str(path))
        return content

    try:
        version = run('uv-version', [uv, '--version'], ROOT).strip()
        if not re.match(r'^uv 0\.11\.16(?:\s|$)', version):
            raise RuntimeError('Expected uv 0.11.16, got ' + version)
        record['uv_version'] = version
        with tempfile.TemporaryDirectory(prefix='sjaracne-uv-install-') as temporary:
            scratch = Path(temporary)
            def fresh_source(destination):
                shutil.copytree(ROOT, destination, ignore=shutil.ignore_patterns(
                    '.git', '.venv', 'venv', '__pycache__', '*.pyc', '*.egg-info',
                    'build', 'dist', '.pytest_cache', '.uv-cache', '.cache',
                    'results', 'local_runs', '*.exe'))
                if any((destination / 'SJARACNe/bin').glob('*.exe')):
                    raise RuntimeError('Fresh source unexpectedly contains native executables')

            source = scratch / 'source'
            fresh_source(source)
            native = source / 'SJARACNe/bin/sjaracne.exe'

            def object_times(source_path=source):
                paths = sorted((source_path / 'SJARACNe/build').rglob('*.o'))
                if not paths:
                    raise RuntimeError('Expected retained native objects from the local source build')
                return {str(path.relative_to(source_path)): path.stat().st_mtime_ns for path in paths}

            def sync(name, venv_path, offline=False, editable=False, source_path=source,
                     path_override=None, expect_failure=False):
                env = environment.copy()
                env['UV_PROJECT_ENVIRONMENT'] = str(venv_path)
                if path_override is not None:
                    env['PATH'] = str(path_override)
                command = [uv, 'sync', '--locked', '--python', python,
                           '--project', source_path, '--cache-dir', cache, '--verbose']
                if not editable:
                    command.append('--no-editable')
                if offline:
                    command.append('--offline')
                return run(name, command, scratch, env, expect_failure=expect_failure)

            def probe(name, venv_path, editable=False, source_path=source):
                executable = venv_path / 'bin/python'
                env = environment.copy()
                env['PATH'] = str(venv_path / 'bin')
                run(name, [executable, '-I', '-c', RUNTIME_PROBE,
                           output / (name + '.json'),
                           str(source_path / 'SJARACNe') if editable else 'installed'], scratch, env)
                observed = json.loads((output / (name + '.json')).read_text(encoding='utf-8'))
                for filename, digest in observed['native_binaries_sha256'].items():
                    if digest != sha256(source_path / 'SJARACNe/bin' / filename):
                        raise RuntimeError('Installed native program differs from rebuilt source: ' + filename)
                return observed

            initial_env = scratch / 'initial-env'
            initial_log = sync('initial-locked-sync', initial_env)
            if not package_was_built(initial_log):
                raise RuntimeError('Fresh source copy did not report building SJARACNe')
            initial = probe('initial-runtime', initial_env)
            initial_times = object_times()
            native_mtime = native.stat().st_mtime_ns
            unchanged_log = sync('unchanged-locked-sync', initial_env)
            if package_was_built(unchanged_log) or object_times() != initial_times \
                    or native.stat().st_mtime_ns != native_mtime:
                raise RuntimeError('Unchanged source unexpectedly rebuilt the package or native objects')

            offline_env = scratch / 'offline-env'
            offline_log = sync('offline-fresh-environment', offline_env, offline=True)
            if package_was_built(offline_log):
                raise RuntimeError('Warm-cache offline install unexpectedly rebuilt unchanged source')
            offline = probe('offline-runtime', offline_env)
            if offline['native_sha256'] != initial['native_sha256']:
                raise RuntimeError('Offline cached install changed the native program')

            # Change only the disposable header, preserving the checkout.
            header = source / 'SJARACNe/src/param.h'
            with header.open('a', encoding='utf-8') as handle:
                handle.write('\n// uv integration probe: native header cache invalidation.\n')
            before_header = object_times()
            header_log = sync('offline-header-rebuild', offline_env, offline=True)
            after_header = object_times()
            if not package_was_built(header_log):
                raise RuntimeError('Native header edit did not invalidate the package cache')
            if set(after_header) != set(before_header) or any(
                    after_header[name] <= before_header[name] for name in before_header):
                raise RuntimeError('Header edit did not recompile every dependent native object')
            header_runtime = probe('header-runtime', offline_env)

            # Installed scientific configuration changes must also invalidate
            # the packaged asset cache, without gratuitous C++ compilation.
            # Use the packaged regulator list, not the hash-protected frozen
            # PCHIP model: even JSON whitespace changes rightly fail its
            # scientific provenance check and are not valid runtime fixtures.
            config = source / 'SJARACNe/config/TF_list.txt'
            with config.open('a', encoding='utf-8') as handle:
                handle.write('\n')
            config_log = sync('offline-config-rebuild', offline_env, offline=True)
            if not package_was_built(config_log) or object_times() != after_header:
                raise RuntimeError('Configuration edit did not rebuild assets without recompiling C++')
            updated = probe('config-runtime', offline_env)
            if updated['config_sha256'] != sha256(config) \
                    or updated['config_sha256'] == initial['config_sha256']:
                raise RuntimeError('Reinstalled configuration does not match the changed source asset')
            # Prove that fresh editable installation fails early with an
            # actionable message for either common missing system prerequisite.
            # Restrict only disposable subprocess PATHs; never alter host tools.
            for missing_tool, available_tool in (('make', 'g++'), ('g++', 'make')):
                label = 'make' if missing_tool == 'make' else 'compiler'
                tool_directory = scratch / ('missing-' + label + '-bin')
                tool_directory.mkdir()
                available = shutil.which(available_tool, path=environment.get('PATH'))
                if not available:
                    raise RuntimeError('Host prerequisite is missing: ' + available_tool)
                (tool_directory / available_tool).symlink_to(available)
                missing_source = scratch / ('missing-' + label + '-source')
                fresh_source(missing_source)
                failure_log = sync('editable-missing-' + label,
                                   scratch / ('missing-' + label + '-env'),
                                   offline=True, editable=True, source_path=missing_source,
                                   path_override=tool_directory, expect_failure=True)
                if ('missing ' + missing_tool + ' on PATH.') not in failure_log \
                        or 'sudo apt install build-essential' not in failure_log:
                    raise RuntimeError('Missing-tool failure lacks the actionable prerequisite message')

            # The default user-facing `uv sync` is editable. It must compile
            # from a truly clean source copy without an earlier Make command.
            editable_source = scratch / 'editable-source'
            fresh_source(editable_source)
            editable_env = scratch / 'editable-env'
            editable_log = sync('editable-initial-locked-sync', editable_env,
                                editable=True, source_path=editable_source)
            if not package_was_built(editable_log):
                raise RuntimeError('Fresh default editable sync did not build SJARACNe')
            editable_initial = probe('editable-initial-runtime', editable_env,
                                     editable=True, source_path=editable_source)
            editable_times = object_times(editable_source)
            unchanged_editable_log = sync('editable-unchanged-locked-sync', editable_env,
                                          offline=True, editable=True,
                                          source_path=editable_source)
            if package_was_built(unchanged_editable_log) \
                    or object_times(editable_source) != editable_times:
                raise RuntimeError('Unchanged editable sync unnecessarily rebuilt the package or objects')

            editable_header = editable_source / 'SJARACNe/src/param.h'
            with editable_header.open('a', encoding='utf-8') as handle:
                handle.write('\n// uv editable probe: native header cache invalidation.\n')
            editable_header_log = sync('editable-offline-header-rebuild', editable_env,
                                      offline=True, editable=True,
                                      source_path=editable_source)
            editable_after_header = object_times(editable_source)
            if not package_was_built(editable_header_log) \
                    or set(editable_times) != set(editable_after_header) or any(
                        editable_after_header[name] <= editable_times[name]
                        for name in editable_times):
                raise RuntimeError('Editable header edit did not rebuild package and all native objects')
            probe('editable-header-runtime', editable_env,
                  editable=True, source_path=editable_source)

            editable_cpp = editable_source / 'SJARACNe/src/util.cpp'
            with editable_cpp.open('a', encoding='utf-8') as handle:
                handle.write('\n// uv editable probe: native source cache invalidation.\n')
            editable_cpp_log = sync('editable-offline-source-rebuild', editable_env,
                                   offline=True, editable=True,
                                   source_path=editable_source)
            editable_after_cpp = object_times(editable_source)
            changed_objects = {name for name in editable_after_cpp
                               if editable_after_cpp[name] > editable_after_header[name]}
            if not package_was_built(editable_cpp_log) \
                    or changed_objects != {'SJARACNe/build/util.o'}:
                raise RuntimeError('Editable source edit did not rebuild its native object and package')
            probe('editable-source-runtime', editable_env,
                  editable=True, source_path=editable_source)

            record['checks'] = {
                'initial_core_noneditable_install': True,
                'unchanged_source_reuses_cached_package': True,
                'offline_fresh_environment': True,
                'native_header_invalidates_package_and_objects': True,
                'configuration_invalidates_packaged_assets': True,
                'runtime_without_uv_on_path': True,
                'missing_make_fails_with_actionable_hint': True,
                'missing_compiler_fails_with_actionable_hint': True,
                'initial_core_editable_install_without_manual_make': True,
                'editable_unchanged_source_reuses_cached_package': True,
                'editable_header_invalidates_package_and_objects': True,
                'editable_source_invalidates_package_and_its_object': True,
                'initial_object_mtimes_ns': initial_times,
                'header_rebuilt_object_mtimes_ns': after_header,
                'editable_initial_object_mtimes_ns': editable_times,
                'editable_header_rebuilt_object_mtimes_ns': editable_after_header,
                'editable_source_rebuilt_object_mtimes_ns': editable_after_cpp,
                'initial_native_sha256': initial['native_sha256'],
                'header_native_sha256': header_runtime['native_sha256'],
                'editable_initial_native_binaries_sha256': editable_initial['native_binaries_sha256'],
            }
        record['complete'] = True
        print('Automatic editable/wheel builds, offline reuse and native/config rebuild probes passed.',
              flush=True)
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
