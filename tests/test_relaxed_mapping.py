"""Exercise real fixmate/markdup ordering with orphan and cross-contig pairs."""
from pathlib import Path
import subprocess
import sys
import pytest
from ipyrad2.mapper import mapper
from ipyrad2.mapper.map_samples_prep import MapperSamplePlan
from ipyrad2.mapper.map_stats import collect_paired_bam_stats, build_map_stats_payload, render_map_stats_payload_report

BIN=Path(sys.prefix)/'bin'
@pytest.mark.parametrize('dedup', [False, True])
def test_relaxed_mapping_preserves_orphan_flags_and_removes_duplicates(tmp_path,monkeypatch,dedup):
    if not (BIN/'samtools').exists():pytest.skip('samtools required')
    sam=tmp_path/'input.sam'
    lines=['@HD\tVN:1.6\tSO:unsorted','@SQ\tSN:a\tLN:1000','@SQ\tSN:b\tLN:1000']
    def rec(name,flag,chrom,pos,mchrom,mpos):
        return '\t'.join(map(str,[name,flag,chrom,pos,60,'100M' if not flag&4 else '*',mchrom,mpos,0,'A'*100,'I'*100,'NM:i:0']))
    for name in ('orphan','orphan_dup'):
        lines += [rec(name,73,'a',50,'=',50),rec(name,133,'a',50,'=',50)]
    for name in ('cross','cross_dup'):
        lines += [rec(name,65,'a',250,'b',250),rec(name,129,'b',250,'a',250)]
    sam.write_text('\n'.join(lines)+'\n')
    monkeypatch.setattr(mapper,'_bwa_mem_cmd',lambda *a:['cat',str(sam)])
    monkeypatch.setattr(mapper,'materialize_sample_plan',lambda *a,**k:((tmp_path/'R1',tmp_path/'R2'),[]))
    out=tmp_path/'mapped';(out/'tmpdir').mkdir(parents=True)
    plan=MapperSamplePlan(output_name='s',source_names=('s',),source_fastqs=((tmp_path/'R1',tmp_path/'R2'),),is_paired_input=True)
    result=mapper._map_sample('s',plan,tmp_path/'ref',out,1,True,False,dedup,False,relaxed_pairs=True)
    rows=[x.split('\t') for x in subprocess.check_output([str(BIN/'samtools'),'view',str(result.bam_path)],text=True).splitlines()]
    assert len(rows)==(3 if dedup else 6)
    assert all(int(x[1])&1 for x in rows)
    assert sum(bool(int(x[1])&8) for x in rows)==(1 if dedup else 2)
    assert all(not int(x[1])&0xF04 for x in rows)
    assert result.nreads_filtered_before_bam_by_mate_unmapped_or_cross_scaffold==0
    assert result.nreads_written_before_duplicate_removal==6
    assert result.nreads_processed == 8
    assert result.nonprimary_alignment_records_removed == 0
    stats = collect_paired_bam_stats(result)
    payload = build_map_stats_payload({'s': stats}, True, relaxed_pairs=True)
    assert payload['relaxed_pairs'] is True
    assert stats['templates_in_final_bam'] == (2 if dedup else 4)
    assert 'regardless of mate location' in render_map_stats_payload_report(payload)


def test_unmate_and_relaxed_pairs_are_mutually_exclusive(tmp_path):
    from ipyrad2.cli.cli_main import setup_parsers
    from ipyrad2.utils.exceptions import IPyradError
    import inspect
    args = setup_parsers().parse_args(['map', '-d', 's_R1.fq', 's_R2.fq', '-r', 'ref.fa'])
    assert args.relaxed_pairs is False
    signature = inspect.signature(mapper.run_mapper)
    kwargs = {k: None for k,p in signature.parameters.items() if p.default is inspect.Parameter.empty}
    kwargs.update(outdir=tmp_path/'must-not-exist', unmate=True, relaxed_pairs=True)
    with pytest.raises(IPyradError, match='cannot be combined'):
        mapper.run_mapper(**kwargs)
    assert not (tmp_path/'must-not-exist').exists()


def test_map_cli_passes_relaxed_pairs_to_mapper(monkeypatch):
    from types import SimpleNamespace
    from ipyrad2.cli import cli_main
    argv = ['map', '-d', 's_R1.fq', 's_R2.fq', '-r', 'ref.fa', '--relaxed-pairs', '-m']
    args = cli_main.setup_parsers().parse_args(argv)
    seen = {}
    monkeypatch.setattr(cli_main.importlib, 'import_module', lambda *a, **k:
                        SimpleNamespace(run_mapper=lambda **kw: seen.update(kw)))
    monkeypatch.setattr(cli_main.sys, 'argv', ['ipyrad2', *argv])
    cli_main.run_subcommand(args, _exit=False)
    assert seen['relaxed_pairs'] and seen['mark_dups_by_coords'] and not seen['unmate']
