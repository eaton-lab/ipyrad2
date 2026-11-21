import gzip

import pytest

from ipyrad2.denovo.sampling import sampled_fastqs
from ipyrad2.utils.exceptions import IPyradError


def _records(start, stop, mate):
    return b"".join(f"@read{i}/{mate}\nACGT\n+\nIIII\n".encode() for i in range(start, stop))


def _ids(path):
    with gzip.open(path, "rb") as handle:
        return [line.split(b"/")[0] for i, line in enumerate(handle) if i % 4 == 0]


@pytest.mark.parametrize("paired", [False, True])
def test_cap_is_reproducible_and_keeps_matched_pairs(tmp_path, paired):
    first = tmp_path / "input_R1.fastq.gz"
    # Multi-member gzip is how pooled RAD observations are staged.
    first.write_bytes(gzip.compress(_records(0, 20, 1)) + gzip.compress(_records(20, 80, 1)))
    second = tmp_path / "input_R2.fastq" if paired else None
    if second:
        second.write_bytes(_records(0, 80, 2))
    scratch = tmp_path / "work"
    scratch.mkdir()
    results = []
    for _ in range(2):
        with sampled_fastqs("sample", first, second, 12, scratch) as (r1, r2, counts):
            ids = _ids(r1)
            assert len(ids) == len(set(ids)) == 12
            assert ids != [f"@read{i}".encode() for i in range(12)]
            assert counts == {"input": 80, "retained": 12}
            if paired:
                assert ids == _ids(r2)
            results.append(ids)
        assert not list(scratch.iterdir())
    assert results[0] == results[1]
    with sampled_fastqs("sample", first, second, 80, scratch) as (r1, r2, counts):
        assert (r1, r2, counts["retained"]) == (first, second, 80)
        assert not list(scratch.iterdir())


def test_failed_sampling_and_clustering_clean_scratch(tmp_path):
    first = tmp_path / "R1.fastq"
    second = tmp_path / "R2.fastq"
    first.write_bytes(_records(0, 20, 1))
    second.write_bytes(_records(1, 21, 2))
    scratch = tmp_path / "work"
    scratch.mkdir()
    with pytest.raises(IPyradError, match="mate IDs differ"):
        with sampled_fastqs("sample", first, second, 12, scratch):
            pytest.fail("mismatched mates reached clustering")
    assert not list(scratch.iterdir())
    second.write_bytes(_records(0, 20, 2))
    with pytest.raises(RuntimeError, match="clustering failed"):
        with sampled_fastqs("sample", first, second, 12, scratch):
            raise RuntimeError("clustering failed")
    assert not list(scratch.iterdir())
