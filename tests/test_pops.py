from __future__ import annotations

import pytest

from ipyrad2.utils.exceptions import IPyradError
from ipyrad2.utils.pops import expand_imap_patterns


def test_expand_imap_patterns_supports_exact_and_glob_matches() -> None:
    expanded, unmatched = expand_imap_patterns(
        {"pop1": ["barbeyi*", "barbeyi-01"], "pop2": ["geyeri-01"]},
        ["barbeyi-01", "barbeyi-02", "geyeri-01"],
    )

    assert expanded == {
        "pop1": ["barbeyi-01", "barbeyi-02"],
        "pop2": ["geyeri-01"],
    }
    assert unmatched == []


def test_expand_imap_patterns_rejects_cross_group_overlap() -> None:
    with pytest.raises(IPyradError, match="assigns sample\\(s\\) multiple times: barbeyi-01"):
        expand_imap_patterns(
            {"pop1": ["barbeyi*"], "pop2": ["barbeyi-01"]},
            ["barbeyi-01", "barbeyi-02"],
        )


def test_expand_imap_patterns_reports_unmatched_entries_when_not_strict() -> None:
    expanded, unmatched = expand_imap_patterns(
        {"pop1": ["barbeyi*", "califo*"], "pop2": ["geyeri-01"]},
        ["barbeyi-01", "geyeri-01"],
        strict_unmatched=False,
    )

    assert expanded == {
        "pop1": ["barbeyi-01"],
        "pop2": ["geyeri-01"],
    }
    assert unmatched == ["califo*"]


def test_expand_imap_patterns_rejects_unmatched_entries_in_strict_mode() -> None:
    with pytest.raises(
        IPyradError,
        match="contains sample names or glob patterns that were not found in the available samples: califo\\*",
    ):
        expand_imap_patterns(
            {"pop1": ["califo*"]},
            ["barbeyi-01"],
        )
