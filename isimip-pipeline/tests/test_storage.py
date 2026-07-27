"""Tests for isimip_pipeline.storage.

These cover key construction, version identifiers and cache mapping — the parts
that must be right for the S3 layout contract in STORAGE.md to hold. They make
no network calls; S3 transfer, the publication gate and the raw-cleanup guards
are exercised separately against the live bucket.
"""

from datetime import datetime, timezone
from pathlib import Path
from unittest import mock

import pytest

from isimip_pipeline import storage


# --- Key construction ---------------------------------------------------


def test_layer_prefix_tcfd():
    assert (
        storage.layer_prefix("wildfire_burntarea_annual")
        == "TCFD/tcfd/layers/wildfire_burntarea_annual"
    )


def test_layer_prefix_water_index_uses_variables_tree():
    """The two products must never share a prefix (CLAUDE.md)."""
    assert (
        storage.layer_prefix("tws", product=storage.PRODUCT_WATER)
        == "TCFD/water-index/variables/tws"
    )


def test_version_prefix():
    assert storage.version_prefix("drought_led_annual", "v2026-07-27_abc1234") == (
        "TCFD/tcfd/layers/drought_led_annual/v2026-07-27_abc1234"
    )


def test_data_key_defaults_to_current():
    assert storage.data_key("cyclone_let_annual", "let_rcp26_processed.nc") == (
        "TCFD/tcfd/layers/cyclone_let_annual/current/let_rcp26_processed.nc"
    )


def test_data_key_versioned_lands_under_data():
    assert storage.data_key(
        "cyclone_let_annual", "let_rcp26_processed.nc", version="v2026-07-27_abc1234"
    ) == (
        "TCFD/tcfd/layers/cyclone_let_annual/v2026-07-27_abc1234/data/"
        "let_rcp26_processed.nc"
    )


def test_raw_prefix_matches_layer_id():
    """Raw and processed prefixes share the layer id — the drift migration fixed."""
    for layer in ("cyclone_let_annual", "fish-tcb_tcb_monthly"):
        assert storage.raw_prefix(layer).endswith(layer)
        assert storage.layer_prefix(layer).endswith(layer)


def test_export_and_reference_keys():
    assert storage.export_prefix("acme", "2026-07-27") == "TCFD/exports/acme/2026-07-27"
    assert storage.reference_key("land_mask.nc") == "TCFD/reference/land_mask.nc"


def test_uri_roundtrips_through_scheme_stripper():
    key = storage.data_key("drought_led_annual", "led_rcp26_processed.nc")
    assert storage._p(storage.uri(key)) == f"{storage.BUCKET}/{key}"


# --- Version identifiers ------------------------------------------------


def test_new_version_format():
    when = datetime(2026, 7, 27, tzinfo=timezone.utc)
    assert storage.new_version(when=when, sha="3412446", dirty=False) == (
        "v2026-07-27_3412446"
    )


def test_new_version_marks_dirty_tree():
    """A dirty tree means the SHA does not identify the code that ran."""
    when = datetime(2026, 7, 27, tzinfo=timezone.utc)
    assert storage.new_version(when=when, sha="3412446", dirty=True) == (
        "v2026-07-27_3412446-dirty"
    )


def test_versions_sort_chronologically():
    a = storage.new_version(when=datetime(2026, 7, 20, tzinfo=timezone.utc),
                            sha="aaaaaaa", dirty=False)
    b = storage.new_version(when=datetime(2026, 7, 27, tzinfo=timezone.utc),
                            sha="bbbbbbb", dirty=False)
    assert sorted([b, a]) == [a, b]


def test_next_available_version_returns_base_when_free():
    with mock.patch.object(storage, "version_exists", return_value=False):
        assert storage.next_available_version("x", base="v2026-07-27_abc") == (
            "v2026-07-27_abc"
        )


def test_next_available_version_bumps_past_taken_ids():
    """Same-day, same-commit reruns must not collide onto a published version."""
    taken = {"v2026-07-27_abc", "v2026-07-27_abc-b"}

    with mock.patch.object(
        storage, "version_exists", side_effect=lambda l, v, p=None: v in taken
    ):
        assert storage.next_available_version("x", base="v2026-07-27_abc") == (
            "v2026-07-27_abc-c"
        )


# --- Cache mapping ------------------------------------------------------


def test_cache_path_mirrors_key_minus_root(monkeypatch):
    monkeypatch.setenv(storage._CACHE_ROOT_ENV, "/tmp/cache-under-test")
    key = "TCFD/tcfd/layers/wildfire_burntarea_annual/v1/data/burntarea_rcp26.nc"
    assert storage.cache_path(key) == Path(
        "/tmp/cache-under-test/tcfd/layers/wildfire_burntarea_annual/v1/data/"
        "burntarea_rcp26.nc"
    )


def test_cache_root_env_override(monkeypatch):
    monkeypatch.setenv(storage._CACHE_ROOT_ENV, "/tmp/elsewhere")
    assert storage.cache_root() == Path("/tmp/elsewhere")


def test_cache_root_defaults_off_repo(monkeypatch):
    """The house rule keeps data off the local/EFS filesystem — cache lives in /tmp."""
    monkeypatch.delenv(storage._CACHE_ROOT_ENV, raising=False)
    assert str(storage.cache_root()).startswith("/tmp")


def test_staging_dir_creates_version_layout(monkeypatch, tmp_path):
    monkeypatch.setenv(storage._CACHE_ROOT_ENV, str(tmp_path))
    stage = storage.staging_dir("some_layer_annual")
    assert sorted(p.name for p in stage.iterdir()) == ["data", "maps", "qa"]


def test_staging_dir_clean_removes_stale_files(monkeypatch, tmp_path):
    """A rerun must not publish leftovers from a previous failed attempt."""
    monkeypatch.setenv(storage._CACHE_ROOT_ENV, str(tmp_path))
    stage = storage.staging_dir("some_layer_annual")
    (stage / "data" / "stale.nc").write_text("stale")

    stage = storage.staging_dir("some_layer_annual", clean=True)
    assert not (stage / "data" / "stale.nc").exists()


# --- Credentials --------------------------------------------------------


def test_use_autorefresh_creds_pops_static_vars(monkeypatch):
    """Pinned static creds cannot refresh and kill long jobs; they must be dropped."""
    for var in ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN"):
        monkeypatch.setenv(var, "pinned")
    monkeypatch.setenv("AWS_CONTAINER_CREDENTIALS_RELATIVE_URI", "/v2/creds")

    storage.use_autorefresh_creds()

    import os

    assert "AWS_ACCESS_KEY_ID" not in os.environ
    assert "AWS_SECRET_ACCESS_KEY" not in os.environ
    assert "AWS_SESSION_TOKEN" not in os.environ
    # The provider endpoint must survive, or creds cannot be re-fetched at all.
    assert os.environ["AWS_CONTAINER_CREDENTIALS_RELATIVE_URI"] == "/v2/creds"


def test_clean_env_strips_static_creds_only(monkeypatch):
    monkeypatch.setenv("AWS_ACCESS_KEY_ID", "pinned")
    monkeypatch.setenv("AWS_DEFAULT_REGION", "us-east-2")

    env = storage.clean_env()

    assert "AWS_ACCESS_KEY_ID" not in env
    assert env["AWS_DEFAULT_REGION"] == "us-east-2"


# --- Guard rails --------------------------------------------------------


def test_publish_derived_rejects_unknown_kind():
    with pytest.raises(ValueError, match="qa.*maps"):
        storage.publish_derived("x", "data", Path("/tmp"))


def test_registry_key_is_control_object():
    """Control/metadata objects are underscore-prefixed, per house convention."""
    assert storage.REGISTRY_KEY.startswith(f"{storage.ROOT}/_")
    assert storage.VERSION_POINTER.startswith("_")
    assert storage.COMPLETE_MARKER.startswith("_")
