from __future__ import annotations

from types import SimpleNamespace

import ipyrad2.utils.profiling as profiling_mod


class _LoggerStub:
    def __init__(self) -> None:
        self.calls: list[tuple[object, ...]] = []

    def debug(self, *args) -> None:
        self.calls.append(args)


def test_profile_stage_logs_elapsed_and_rss_metrics(monkeypatch) -> None:
    logger = _LoggerStub()
    times = iter([10.0, 12.5])
    rss_values = iter([10 * 1024 * 1024, 16 * 1024 * 1024])
    usages = iter(
        [
            SimpleNamespace(ru_maxrss=1024),
            SimpleNamespace(ru_maxrss=2048),
            SimpleNamespace(ru_maxrss=3072),
            SimpleNamespace(ru_maxrss=4096),
        ]
    )

    monkeypatch.setattr(profiling_mod, "logger", logger)
    monkeypatch.setattr(profiling_mod.sys, "platform", "linux")
    monkeypatch.setattr(profiling_mod.time, "perf_counter", lambda: next(times))
    monkeypatch.setattr(
        profiling_mod,
        "_current_self_rss_bytes",
        lambda: next(rss_values),
    )
    monkeypatch.setattr(
        profiling_mod.resource,
        "getrusage",
        lambda who: next(usages),
    )

    with profiling_mod.profile_stage("variant calling"):
        pass

    assert logger.calls == [
        (
            "stage profile {}: elapsed={:.2f}s, self_rss_start={}, self_rss_end={}, self_rss_delta={}, self_peak_delta={}, child_peak_delta={}, self_peak_total={}, child_peak_total={}",
            "variant calling",
            2.5,
            "10.0 MiB",
            "16.0 MiB",
            "6.0 MiB",
            "2.0 MiB",
            "2.0 MiB",
            "3.0 MiB",
            "4.0 MiB",
        )
    ]


def test_profile_stage_logs_na_for_unavailable_current_rss(monkeypatch) -> None:
    logger = _LoggerStub()
    times = iter([1.0, 2.0])
    usages = iter(
        [
            SimpleNamespace(ru_maxrss=1024),
            SimpleNamespace(ru_maxrss=1024),
            SimpleNamespace(ru_maxrss=2048),
            SimpleNamespace(ru_maxrss=3072),
        ]
    )

    monkeypatch.setattr(profiling_mod, "logger", logger)
    monkeypatch.setattr(profiling_mod.sys, "platform", "linux")
    monkeypatch.setattr(profiling_mod.time, "perf_counter", lambda: next(times))
    monkeypatch.setattr(
        profiling_mod,
        "_current_self_rss_bytes",
        lambda: None,
    )
    monkeypatch.setattr(
        profiling_mod.resource,
        "getrusage",
        lambda who: next(usages),
    )

    with profiling_mod.profile_stage("sample mask building"):
        pass

    assert logger.calls == [
        (
            "stage profile {}: elapsed={:.2f}s, self_rss_start={}, self_rss_end={}, self_rss_delta={}, self_peak_delta={}, child_peak_delta={}, self_peak_total={}, child_peak_total={}",
            "sample mask building",
            1.0,
            "n/a",
            "n/a",
            "n/a",
            "1.0 MiB",
            "2.0 MiB",
            "2.0 MiB",
            "3.0 MiB",
        )
    ]


def test_format_peak_rss_uses_platform_specific_units(monkeypatch) -> None:
    monkeypatch.setattr(profiling_mod.sys, "platform", "linux")
    assert profiling_mod._format_peak_rss(1024) == "1.0 MiB"

    monkeypatch.setattr(profiling_mod.sys, "platform", "darwin")
    assert profiling_mod._format_peak_rss(1024 * 1024) == "1.0 MiB"
