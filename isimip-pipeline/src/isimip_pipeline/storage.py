"""S3-canonical storage for the TCFD pipeline.

S3 holds the data; the local filesystem is ephemeral cache only. This module is
the **single place S3 keys are constructed** — no other module or script should
build a key by hand. See STORAGE.md for the layout contract.

Layout (under ``s3://climate-ai-data-science-shiny-app-data/TCFD/``)::

    tcfd/layers/{layer_id}/{version}/{data,qa,maps}/  layer.json  _COMPLETE.json
    tcfd/layers/{layer_id}/current/                   # copy of active version's data/
    tcfd/layers/{layer_id}/_VERSION.json
    water-index/variables/{variable}/{version}/…
    exports/{customer}/{run_date}/
    reference/
    raw/isimip/{layer_id}/
    _registry/layers.json

Credentials: this module deliberately drops static ``AWS_*`` env vars so botocore
uses the auto-refreshing SageMaker container-credential provider. Pinning a static
copy of a ~1h task token is the org's costliest documented failure mode (long jobs
die mid-run with ``OSError: Bad Request``). See STORAGE.md § Credentials.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

import s3fs

# --- Location -------------------------------------------------------------

BUCKET = "climate-ai-data-science-shiny-app-data"
ROOT = "TCFD"
REGION = "us-east-2"

PRODUCT_TCFD = "tcfd"
PRODUCT_WATER = "water-index"

COMPLETE_MARKER = "_COMPLETE.json"
VERSION_POINTER = "_VERSION.json"
REGISTRY_KEY = f"{ROOT}/_registry/layers.json"

#: Subdirectories inside a published version prefix.
VERSION_KINDS = ("data", "qa", "maps")

_CACHE_ROOT_ENV = "TCFD_CACHE_ROOT"
_DEFAULT_CACHE_ROOT = "/tmp/tcfd-cache"

_STATIC_CRED_VARS = (
    "AWS_ACCESS_KEY_ID",
    "AWS_SECRET_ACCESS_KEY",
    "AWS_SESSION_TOKEN",
    "AWS_SECURITY_TOKEN",
)


# --- Credentials ----------------------------------------------------------


def use_autorefresh_creds() -> None:
    """Drop static AWS creds so botocore uses the auto-refreshing provider.

    SageMaker task-role tokens are temporary (~1 hour). A static copy pinned into
    ``AWS_*`` env vars cannot refresh, so any S3 call after expiry fails — the
    documented failure mode is a long job dying ~90 minutes in. Popping the static
    vars (``AWS_CONTAINER_CREDENTIALS_RELATIVE_URI`` stays) lets the container
    provider re-fetch between requests.

    Idempotent. Re-assert before each fresh round of S3 I/O in a long-running job.
    """
    for var in _STATIC_CRED_VARS:
        os.environ.pop(var, None)
    s3fs.S3FileSystem.clear_instance_cache()


def clean_env() -> Dict[str, str]:
    """Return ``os.environ`` stripped of ``AWS_*`` static creds.

    For ``aws s3 cp`` subprocesses, which re-resolve IAM per invocation. Every
    such subprocess — read or write — must use this; a silent empty-result
    corruption in a sibling project came from omitting it once.
    """
    return {k: v for k, v in os.environ.items() if k not in _STATIC_CRED_VARS}


_fs: Optional[s3fs.S3FileSystem] = None


def s3_filesystem(refresh: bool = False) -> s3fs.S3FileSystem:
    """Return a cached :class:`s3fs.S3FileSystem` using auto-refreshing IAM creds.

    Args:
        refresh: Re-assert credential handling and rebuild the filesystem. Use
            this between long-running phases of a job.
    """
    global _fs
    if _fs is None or refresh:
        use_autorefresh_creds()
        _fs = s3fs.S3FileSystem(anon=False)
    return _fs


def _p(uri: str) -> str:
    """Normalize a URI to the scheme-less ``bucket/key`` form s3fs expects."""
    return uri.replace("s3://", "").lstrip("/")


def _guard_write_key(key: str) -> str:
    """Reject any write destination outside ``{ROOT}/``.

    Writes are confined to ``s3://{BUCKET}/{ROOT}/`` by policy: every other prefix
    in this bucket, and every other bucket, belongs to a different live product and
    is read-only. ``BUCKET`` is already pinned as a module constant, but the ``key``
    passed to :func:`push` / :func:`push_dir` / :func:`write_json` is
    caller-supplied, so a malformed key could otherwise land in a sibling prefix.

    Args:
        key: Bucket-relative destination key.

    Returns:
        The key, stripped of any leading slash.

    Raises:
        ValueError: If the key does not sit under ``{ROOT}/``, or escapes it via
            ``..`` path traversal.
    """
    clean = str(key).lstrip("/")
    if not clean:
        raise ValueError("refusing to write to an empty key")
    if ".." in Path(clean).parts:
        raise ValueError(f"refusing to write to a key containing '..': {key!r}")
    if clean != ROOT and not clean.startswith(f"{ROOT}/"):
        raise ValueError(
            f"refusing to write outside s3://{BUCKET}/{ROOT}/ — got {key!r}. "
            "Build every key through this module (data_key, raw_prefix, "
            "raw_inhouse_prefix, reference_key, …); never hand-assemble one."
        )
    return clean


def uri(key: str) -> str:
    """Return the full ``s3://`` URI for a bucket-relative key."""
    return f"s3://{BUCKET}/{key.lstrip('/')}"


def exists(key: str) -> bool:
    """Return True if an object or prefix exists."""
    return s3_filesystem().exists(_p(f"{BUCKET}/{key}"))


# --- Versions -------------------------------------------------------------


def git_sha(short: int = 7) -> str:
    """Return the current short git SHA, or ``"nogit"`` outside a repository."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", f"--short={short}", "HEAD"],
            cwd=Path(__file__).resolve().parents[3],
            capture_output=True,
            text=True,
            check=True,
        )
        return out.stdout.strip() or "nogit"
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return "nogit"


def git_dirty() -> bool:
    """Return True if the working tree has uncommitted changes."""
    try:
        out = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=Path(__file__).resolve().parents[3],
            capture_output=True,
            text=True,
            check=True,
        )
        return bool(out.stdout.strip())
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return False


def new_version(
    when: Optional[datetime] = None,
    sha: Optional[str] = None,
    dirty: Optional[bool] = None,
) -> str:
    """Build a version identifier: ``v{YYYY-MM-DD}_{short-git-sha}``.

    Chronologically sortable and commit-traceable; same-day republishes from
    different commits never collide. A dirty working tree appends ``-dirty``,
    because the SHA alone does not then identify the code that ran.

    Same-day republishes from the *same* commit DO collide by construction — use
    :func:`next_available_version` to disambiguate, or let
    :func:`publish_layer_version` refuse the overwrite.

    Args:
        when: UTC timestamp to date-stamp with. Defaults to now.
        sha: Short git SHA. Defaults to the current HEAD.
        dirty: Override dirty-tree detection.
    """
    when = when or datetime.now(timezone.utc)
    suffix = "-dirty" if (git_dirty() if dirty is None else dirty) else ""
    return f"v{when.strftime('%Y-%m-%d')}_{sha or git_sha()}{suffix}"


def version_exists(layer_id: str, version: str, product: str = PRODUCT_TCFD) -> bool:
    """Return True if a version prefix already holds a completion marker."""
    return exists(f"{version_prefix(layer_id, version, product)}/{COMPLETE_MARKER}")


def next_available_version(
    layer_id: str,
    product: str = PRODUCT_TCFD,
    base: Optional[str] = None,
) -> str:
    """Return ``base`` if unpublished, else the first free ``-b``, ``-c``, … variant.

    Lets a same-day, same-commit rerun publish alongside its predecessor rather
    than overwriting it.
    """
    base = base or new_version()
    if not version_exists(layer_id, base, product):
        return base
    for letter in "bcdefghijklmnopqrstuvwxyz":
        candidate = f"{base}-{letter}"
        if not version_exists(layer_id, candidate, product):
            return candidate
    raise RuntimeError(f"{layer_id}: exhausted version suffixes for {base}")


# --- Key construction ----------------------------------------------------


def layer_prefix(layer_id: str, product: str = PRODUCT_TCFD) -> str:
    """Prefix holding every version of one layer."""
    if product == PRODUCT_WATER:
        return f"{ROOT}/{PRODUCT_WATER}/variables/{layer_id}"
    return f"{ROOT}/{product}/layers/{layer_id}"


def version_prefix(layer_id: str, version: str, product: str = PRODUCT_TCFD) -> str:
    """Prefix of one immutable published version."""
    return f"{layer_prefix(layer_id, product)}/{version}"


def current_prefix(layer_id: str, product: str = PRODUCT_TCFD) -> str:
    """Stable read path for the active version's data files."""
    return f"{layer_prefix(layer_id, product)}/current"


def version_pointer_key(layer_id: str, product: str = PRODUCT_TCFD) -> str:
    """Key of the per-layer ``_VERSION.json`` pointer."""
    return f"{layer_prefix(layer_id, product)}/{VERSION_POINTER}"


def data_key(
    layer_id: str,
    filename: str,
    version: str = "current",
    product: str = PRODUCT_TCFD,
) -> str:
    """Key of a single data file, in ``current/`` or inside a specific version.

    Args:
        layer_id: Canonical layer id, e.g. ``wildfire_burntarea_annual``.
        filename: Bare filename, e.g. ``burntarea_rcp26_processed.nc``.
        version: ``"current"`` (default) or a version id from :func:`new_version`.
        product: ``"tcfd"`` or ``"water-index"``.
    """
    if version == "current":
        return f"{current_prefix(layer_id, product)}/{filename}"
    return f"{version_prefix(layer_id, version, product)}/data/{filename}"


def raw_prefix(layer_id: str) -> str:
    """Staging prefix for downloaded ISIMIP members of one layer."""
    return f"{ROOT}/raw/isimip/{layer_id}"


def raw_inhouse_prefix(dataset_id: str) -> str:
    """Staging prefix for a raw member set built in-house rather than pulled from ISIMIP.

    Sits alongside ``raw/isimip/`` so the provenance is visible in the key itself:
    these members are derived from internal archives (e.g. the bias-corrected CMIP6
    ``ws_max`` store), not downloaded from the ISIMIP file server, so they are not
    re-fetchable by URL and ``cleanup_raw``'s source-URL safety check does not apply.

    Args:
        dataset_id: Canonical dataset id, e.g. ``windydays_wsmax17p5_annual``.
    """
    return f"{ROOT}/raw/in-house/{dataset_id}"


def export_prefix(customer: str, run_date: str) -> str:
    """Prefix for one customer's Export-Key CSV delivery."""
    return f"{ROOT}/exports/{customer}/{run_date}"


def reference_key(filename: str) -> str:
    """Key of a shared reference input (land mask, region polygons, …)."""
    return f"{ROOT}/reference/{filename}"


def dev_prefix(name: str) -> str:
    """Throwaway prefix for smoke tests and scratch work — safe to delete."""
    return f"{ROOT}/_dev/{name}"


# --- Local cache ---------------------------------------------------------


def cache_root() -> Path:
    """Root of the ephemeral local cache (``$TCFD_CACHE_ROOT`` or ``/tmp/tcfd-cache``)."""
    return Path(os.environ.get(_CACHE_ROOT_ENV, _DEFAULT_CACHE_ROOT))


def cache_path(key: str) -> Path:
    """Local cache path mirroring a bucket-relative key.

    The cache mirrors the S3 key path exactly, minus the ``TCFD/`` root, so the
    cache-to-S3 mapping is trivially reversible.
    """
    rel = key.lstrip("/")
    if rel.startswith(f"{ROOT}/"):
        rel = rel[len(ROOT) + 1 :]
    return cache_root() / rel


def clear_cache() -> None:
    """Delete the entire local cache. Never loses data — S3 is canonical."""
    root = cache_root()
    if root.exists():
        shutil.rmtree(root)


def staging_dir(layer_id: str, clean: bool = True) -> Path:
    """Return a local directory laid out like a version prefix, for building one.

    Creates ``data/``, ``qa/`` and ``maps/`` subdirectories. Processors write
    their outputs here, then hand the directory to
    :func:`publish_layer_version`.

    Args:
        layer_id: Canonical layer id.
        clean: Remove any previous staging contents first, so a rerun cannot
            publish stale files left over from an earlier attempt.
    """
    stage = cache_root() / "_staging" / layer_id
    if clean and stage.exists():
        shutil.rmtree(stage)
    for kind in VERSION_KINDS:
        (stage / kind).mkdir(parents=True, exist_ok=True)
    return stage


# --- Checksums ----------------------------------------------------------


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    """Return the hex sha256 of a local file, read in chunks."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


# --- Transfer -----------------------------------------------------------


def pull(key: str, force: bool = False) -> Path:
    """Download one object into the local cache and return its path.

    A cached file is reused when its size matches the S3 object. Pass
    ``force=True`` to re-download unconditionally.

    Args:
        key: Bucket-relative key.
        force: Ignore any cached copy.

    Returns:
        Local path to the cached file.

    Raises:
        FileNotFoundError: If the object does not exist.
    """
    fs = s3_filesystem()
    src = _p(f"{BUCKET}/{key}")
    if not fs.exists(src):
        raise FileNotFoundError(f"s3://{src}")

    dest = cache_path(key)
    if dest.exists() and not force:
        if dest.stat().st_size == fs.info(src)["size"]:
            return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(f"{dest.name}.downloading-{uuid.uuid4().hex[:8]}")
    try:
        fs.get(src, str(tmp))
        tmp.replace(dest)
    finally:
        if tmp.exists():
            tmp.unlink()
    return dest


def push(local: Path, key: str, atomic: bool = True) -> Dict[str, Any]:
    """Upload one file and return its manifest entry.

    Args:
        local: Local file to upload.
        key: Bucket-relative destination key.
        atomic: Upload to a temp key then server-side move into place, so readers
            never observe a partial object. Disable for large raw staging files,
            where the extra copy is wasteful and no completeness gate depends on
            the object.

    Returns:
        ``{"key", "bytes", "sha256"}`` — ready to drop into a manifest.
    """
    local = Path(local)
    key = _guard_write_key(key)
    fs = s3_filesystem()
    dest = _p(f"{BUCKET}/{key}")

    if atomic:
        tmp = f"{dest}.writing-{uuid.uuid4().hex[:12]}"
        fs.put(str(local), tmp)
        if fs.exists(dest):
            fs.rm(dest)
        fs.mv(tmp, dest)
    else:
        fs.put(str(local), dest)

    return {"key": key, "bytes": local.stat().st_size, "sha256": sha256_file(local)}


def push_dir(local_dir: Path, prefix: str, atomic: bool = True) -> List[Dict[str, Any]]:
    """Upload a directory tree and return manifest entries for every file.

    Manifest ``key`` values are relative to ``prefix``, which is what
    :func:`write_complete_marker` expects.
    """
    local_dir = Path(local_dir)
    _guard_write_key(prefix)  # fail fast, even when the tree is empty
    entries: List[Dict[str, Any]] = []
    for path in sorted(p for p in local_dir.rglob("*") if p.is_file()):
        rel = path.relative_to(local_dir).as_posix()
        entry = push(path, f"{prefix.rstrip('/')}/{rel}", atomic=atomic)
        entry["key"] = rel
        entries.append(entry)
    return entries


def pull_prefix(prefix: str, force: bool = False) -> Path:
    """Download every object under a prefix into the cache; return the local dir."""
    fs = s3_filesystem()
    src = _p(f"{BUCKET}/{prefix}")
    keys = [k for k in fs.find(src) if not k.endswith("/")]
    if not keys:
        raise FileNotFoundError(f"s3://{src} (empty or missing)")
    for k in keys:
        pull(k[len(BUCKET) + 1 :], force=force)
    return cache_path(prefix)


def stage_raw(layer_id: str, pattern: str = "*", force: bool = False) -> List[Path]:
    """Download a layer's raw ISIMIP members into the cache and return their paths.

    NetCDF libraries cannot read lazily from S3, so members are materialized
    locally. Already-cached members are reused unless ``force`` is set.

    Args:
        layer_id: Canonical layer id.
        pattern: Filename glob applied to the raw prefix, e.g. ``"*_2006_2099.nc4"``.
        force: Re-download even if cached.

    Returns:
        Sorted local paths. Empty if the raw prefix holds nothing matching.
    """
    fs = s3_filesystem()
    prefix = _p(f"{BUCKET}/{raw_prefix(layer_id)}")
    if not fs.exists(prefix):
        return []
    matches = [k for k in fs.glob(f"{prefix}/{pattern}") if not k.endswith("/")]
    return sorted(pull(k[len(BUCKET) + 1 :], force=force) for k in matches)


def ingest_raw(
    files: Sequence[Path],
    layer_id: str,
    source_urls: Optional[Dict[str, str]] = None,
) -> List[Dict[str, Any]]:
    """Upload downloaded ISIMIP members to a layer's raw staging prefix.

    Uploads are non-atomic: raw objects are large and no completeness gate
    depends on them, so the extra server-side copy would be wasted I/O.

    Args:
        files: Local member files.
        layer_id: Canonical layer id.
        source_urls: Optional ``{filename: url}`` map. Recording the origin URL is
            what later makes :func:`cleanup_raw` safe, so pass it whenever known.

    Returns:
        Manifest entries (``key``, ``bytes``, ``sha256``, and ``source_url`` when
        supplied), suitable for ``layer.json``'s ``inputs.files``.
    """
    prefix = raw_prefix(layer_id)
    entries: List[Dict[str, Any]] = []
    for path in sorted(Path(p) for p in files):
        entry = push(path, f"{prefix}/{path.name}", atomic=False)
        entry["name"] = path.name
        entry.pop("key", None)
        if source_urls and path.name in source_urls:
            entry["source_url"] = source_urls[path.name]
        entries.append(entry)
    return entries


def publish_derived(layer_id: str, kind: str, local_dir: Path,
                    version: Optional[str] = None,
                    product: str = PRODUCT_TCFD) -> str:
    """Upload regenerable derived artifacts (``qa`` or ``maps``) into a version.

    These live inside the version prefix — so evidence stays pinned to the data
    version it describes — but sit **outside** the ``_COMPLETE.json`` gate, which
    covers only the immutable data contract. That lets QA reports and maps be
    regenerated as visualization code improves without invalidating the gate or
    forcing a new data version.

    Args:
        layer_id: Canonical layer id.
        kind: ``"qa"`` or ``"maps"``.
        local_dir: Directory whose contents to upload.
        version: Target version. Defaults to the layer's current version.
        product: ``"tcfd"`` or ``"water-index"``.

    Returns:
        The prefix written to.
    """
    if kind not in ("qa", "maps"):
        raise ValueError(f"kind must be 'qa' or 'maps', got {kind!r}")
    version = version or resolve_current(layer_id, product)
    prefix = f"{version_prefix(layer_id, version, product)}/{kind}"
    push_dir(Path(local_dir), prefix)
    return prefix


def read_json(key: str) -> Dict[str, Any]:
    """Read and parse a JSON object from S3."""
    fs = s3_filesystem()
    with fs.open(_p(f"{BUCKET}/{key}"), "rb") as f:
        return json.loads(f.read().decode("utf-8"))


def write_json(obj: Dict[str, Any], key: str) -> str:
    """Write a dict to S3 as indented JSON. Returns the key."""
    key = _guard_write_key(key)
    fs = s3_filesystem()
    body = json.dumps(obj, indent=2, default=str).encode("utf-8")
    dest = _p(f"{BUCKET}/{key}")
    tmp = f"{dest}.writing-{uuid.uuid4().hex[:12]}"
    with fs.open(tmp, "wb") as f:
        f.write(body)
    if fs.exists(dest):
        fs.rm(dest)
    fs.mv(tmp, dest)
    return key


# --- Publication gate ---------------------------------------------------


def write_complete_marker(
    prefix: str,
    artifacts: Sequence[Dict[str, Any]],
    extra: Optional[Dict[str, Any]] = None,
) -> str:
    """Write ``_COMPLETE.json`` — the sha256 gate — as the LAST step of a publish.

    Its presence and validity is what makes a version consumable. Anything that
    reads a version must call :func:`verify_complete` first.

    Args:
        prefix: Version prefix the marker lives in.
        artifacts: Manifest entries from :func:`push` / :func:`push_dir`, with
            ``key`` relative to ``prefix``.
        extra: Additional fields to record alongside the artifact table.
    """
    marker = {
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "artifacts": {
            a["key"]: {"sha256": a["sha256"], "bytes": a["bytes"]} for a in artifacts
        },
        **(extra or {}),
    }
    return write_json(marker, f"{prefix.rstrip('/')}/{COMPLETE_MARKER}")


def verify_complete(prefix: str, require: Iterable[str] = ()) -> Dict[str, Any]:
    """Verify a version's completion marker, re-hashing every artifact.

    Args:
        prefix: Version prefix to verify.
        require: Artifact keys that must be present in the marker.

    Returns:
        The parsed marker.

    Raises:
        FileNotFoundError: No marker — the version is in-flight or pre-dates the
            convention. Never treat such a prefix as data.
        ValueError: A required artifact is missing, or a sha256 does not match
            (mixed-version or corrupted publication).
    """
    fs = s3_filesystem()
    prefix = prefix.rstrip("/")
    key = f"{prefix}/{COMPLETE_MARKER}"
    if not exists(key):
        raise FileNotFoundError(
            f"s3://{BUCKET}/{prefix}: no {COMPLETE_MARKER} — publication "
            f"incomplete; refuse to consume."
        )
    marker = read_json(key)
    artifacts = marker.get("artifacts", {})

    for name in require:
        if name not in artifacts:
            raise ValueError(f"s3://{BUCKET}/{prefix}: marker lacks artifact {name!r}")

    for name, info in artifacts.items():
        body = fs.cat(_p(f"{BUCKET}/{prefix}/{name}"))
        got = hashlib.sha256(body).hexdigest()
        if got != info["sha256"]:
            raise ValueError(
                f"s3://{BUCKET}/{prefix}/{name}: sha256 mismatch vs marker — "
                f"mixed-version or corrupted publication"
            )
    return marker


# --- Registry and version pointer --------------------------------------


def read_registry() -> Dict[str, Any]:
    """Read ``_registry/layers.json``, or an empty registry if absent."""
    if not exists(REGISTRY_KEY):
        return {"schema_version": 1, "updated_utc": None, "layers": {}}
    return read_json(REGISTRY_KEY)


def update_registry_entry(layer_id: str, entry: Dict[str, Any]) -> str:
    """Merge one layer's entry into the registry and write it back."""
    reg = read_registry()
    reg.setdefault("layers", {})
    current = reg["layers"].get(layer_id, {})
    current.update(entry)
    reg["layers"][layer_id] = current
    reg["updated_utc"] = datetime.now(timezone.utc).isoformat()
    return write_json(reg, REGISTRY_KEY)


def resolve_current(layer_id: str, product: str = PRODUCT_TCFD) -> str:
    """Return the active version id for a layer.

    Raises:
        FileNotFoundError: If the layer has no published version.
    """
    key = version_pointer_key(layer_id, product)
    if not exists(key):
        raise FileNotFoundError(
            f"{layer_id}: no {VERSION_POINTER} — layer has no published version"
        )
    return read_json(key)["current"]


# --- Publish ------------------------------------------------------------


def publish_layer_version(
    layer_id: str,
    local_dir: Path,
    manifest: Dict[str, Any],
    product: str = PRODUCT_TCFD,
    version: Optional[str] = None,
    make_current: bool = True,
    on_exists: str = "error",
) -> str:
    """Publish a layer version, in the order that keeps consumers safe.

    ``local_dir`` is staged like a version prefix::

        local_dir/data/{var}_{scenario}_processed.nc
        local_dir/qa/…
        local_dir/maps/…

    Steps, per STORAGE.md § Publish protocol:

    1. Upload ``data/``, ``qa/``, ``maps/`` (atomic temp-key then move).
    2. Upload ``layer.json``.
    3. Upload ``_COMPLETE.json`` — marks the version consumable.
    4. Server-side copy ``data/*`` into ``current/``.
    5. Write ``_VERSION.json``.
    6. Update ``_registry/layers.json``.

    A failure before step 3 leaves an orphaned version that consumers ignore.

    Args:
        layer_id: Canonical layer id.
        local_dir: Staged version directory.
        manifest: Science manifest; ``layer_id``, ``product``, ``version``,
            ``created_utc`` and ``git_commit`` are filled in if absent.
        product: ``"tcfd"`` or ``"water-index"``.
        version: Version id. Defaults to :func:`new_version`.
        make_current: Copy ``data/`` into ``current/`` and move the pointer.
        on_exists: What to do when the target version is already published —
            ``"error"`` (default, protects immutability), ``"bump"`` (publish
            alongside as ``-b``, ``-c``, … via :func:`next_available_version`),
            or ``"overwrite"`` (replace it; only for a known-bad publish).

    Returns:
        The published version id.

    Raises:
        FileExistsError: If the version is already published and
            ``on_exists="error"``.
    """
    local_dir = Path(local_dir)
    if not (local_dir / "data").is_dir():
        raise ValueError(f"{local_dir}: expected a data/ subdirectory")

    version = version or new_version()
    if version_exists(layer_id, version, product):
        if on_exists == "error":
            raise FileExistsError(
                f"{layer_id} {version} is already published — versions are "
                f"immutable. Pass on_exists='bump' to publish alongside it, or "
                f"on_exists='overwrite' to replace it."
            )
        if on_exists == "bump":
            version = next_available_version(layer_id, product, base=version)
        elif on_exists != "overwrite":
            raise ValueError(f"on_exists must be error|bump|overwrite, got {on_exists!r}")

    vprefix = version_prefix(layer_id, version, product)
    fs = s3_filesystem()

    # 1. artifacts. Only data/ enters the completeness gate — qa/ and maps/ are
    #    regenerable derived output (see publish_derived).
    entries: List[Dict[str, Any]] = []
    for entry in push_dir(local_dir / "data", f"{vprefix}/data"):
        entry["key"] = f"data/{entry['key']}"
        entries.append(entry)
    if not entries:
        raise ValueError(f"{local_dir}/data is empty — nothing to publish")

    for kind in ("qa", "maps"):
        kind_dir = local_dir / kind
        if kind_dir.is_dir() and any(kind_dir.rglob("*")):
            push_dir(kind_dir, f"{vprefix}/{kind}")

    # 2. science manifest
    manifest = dict(manifest)
    manifest.setdefault("schema_version", 1)
    manifest.setdefault("layer_id", layer_id)
    manifest.setdefault("product", product)
    manifest.setdefault("version", version)
    manifest.setdefault("created_utc", datetime.now(timezone.utc).isoformat())
    manifest.setdefault("git_commit", git_sha())

    manifest_local = cache_path(f"{vprefix}/layer.json")
    manifest_local.parent.mkdir(parents=True, exist_ok=True)
    manifest_local.write_text(json.dumps(manifest, indent=2, default=str))
    entries.append(push(manifest_local, f"{vprefix}/layer.json"))
    entries[-1]["key"] = "layer.json"

    # 3. gate
    write_complete_marker(vprefix, entries, extra={"version": version, "layer_id": layer_id})

    if make_current:
        # 4. current/ copy (server-side)
        cprefix = current_prefix(layer_id, product)
        for entry in entries:
            if not entry["key"].startswith("data/"):
                continue
            name = entry["key"][len("data/") :]
            fs.copy(
                _p(f"{BUCKET}/{vprefix}/data/{name}"),
                _p(f"{BUCKET}/{cprefix}/{name}"),
            )

        # 5. pointer
        history = []
        pointer_key = version_pointer_key(layer_id, product)
        if exists(pointer_key):
            history = read_json(pointer_key).get("history", [])
        if version not in history:
            history.append(version)
        write_json(
            {
                "layer_id": layer_id,
                "current": version,
                "updated_utc": datetime.now(timezone.utc).isoformat(),
                "history": history,
            },
            pointer_key,
        )

        # 6. registry
        update_registry_entry(
            layer_id,
            {
                "product": product,
                "current": version,
                "versions": history,
                "variable": manifest.get("variable"),
                "units": manifest.get("units"),
                "scenarios": manifest.get("scenarios"),
                "percentile_direction": (manifest.get("decisions") or {}).get(
                    "percentile_direction"
                ),
            },
        )

    return version


# --- Raw staging cleanup -----------------------------------------------


def cleanup_raw(
    layer_id: str,
    version: Optional[str] = None,
    product: str = PRODUCT_TCFD,
    dry_run: bool = True,
) -> List[str]:
    """Delete a layer's raw staging prefix, but only once its outputs are safe.

    Refuses to delete unless the published version's ``_COMPLETE.json`` verifies
    and its ``layer.json`` records a ``source_url`` for every input — otherwise
    the raw data could not be re-fetched.

    Args:
        layer_id: Canonical layer id.
        version: Version to validate against. Defaults to the current one.
        product: ``"tcfd"`` or ``"water-index"``.
        dry_run: Report what would be deleted without deleting. Default True.

    Returns:
        Keys deleted, or that would be deleted under ``dry_run``.

    Raises:
        FileNotFoundError, ValueError: If the safety checks do not pass.
    """
    version = version or resolve_current(layer_id, product)
    vprefix = version_prefix(layer_id, version, product)
    verify_complete(vprefix, require=["layer.json"])

    manifest = read_json(f"{vprefix}/layer.json")
    files = (manifest.get("inputs") or {}).get("files") or []
    if not files:
        raise ValueError(
            f"{layer_id} {version}: layer.json records no inputs.files — refusing "
            f"to delete raw that could not be re-fetched"
        )
    missing = [f.get("name", "?") for f in files if not f.get("source_url")]
    if missing:
        raise ValueError(
            f"{layer_id} {version}: {len(missing)} input(s) lack source_url "
            f"(e.g. {missing[0]}) — refusing to delete unrecoverable raw"
        )

    fs = s3_filesystem()
    prefix = _p(f"{BUCKET}/{raw_prefix(layer_id)}")
    keys = [k for k in fs.find(prefix) if not k.endswith("/")]
    if not dry_run and keys:
        fs.rm(keys)
    return [k[len(BUCKET) + 1 :] for k in keys]
